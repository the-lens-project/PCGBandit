import json
import os
import re
import sys
from collections import defaultdict
from glob import glob
import h5py
import numpy as np

# scipy < 1.8 has no sparse *_array classes, so fall back to the *_matrix
# equivalents to run on the testing interpreter (scipy 1.7.x). Only the shared
# interface is used here -- the 3-tuple/(data,(row,col)) constructors, tocsr,
# tocoo, .data/.row/.col, transpose and @ -- so the two are interchangeable.
try:
    from scipy.sparse import csr_array, coo_array
except ImportError:  # scipy < 1.8
    from scipy.sparse import csr_matrix as csr_array, coo_matrix as coo_array


# Quantities stored as integer cell indices (labels).
INT_QUANTS = {'lowerAddr', 'upperAddr', 'interfaceRow', 'interfaceCol'}
# Matrix values (diag/lower/interfaceLower) and addressing (lowerAddr/upperAddr/
# interfaceRow/interfaceCol), deduplicated on the .H side against the most recent
# dump: a file is written only when it changes, so an absent file means "unchanged
# since the previous dump" and the most recently loaded value is carried forward.
# The carry-forward state holds exactly these; b/sol/init are written every solve
# and read straight from the solve directory when a system is assembled.
DEDUP_QUANTS = ['diag', 'lower', 'lowerAddr', 'upperAddr',
                'interfaceRow', 'interfaceCol', 'interfaceLower']


def _parse_block(buf, dtype):
    """Parse a binary OpenFOAM `<count>\\n(<data>)` payload. The data begins at
    the first '(' (neither a FoamFile header nor the count contains one) and the
    count is the integer immediately before it. An empty field is written as a
    bare '0' with no parens, which is the p < 0 case."""
    p = buf.find(b'(')
    if p < 0:
        return np.empty(0, dtype=dtype)
    count = int(re.search(rb'(\d+)\s*$', buf[:p]).group(1))
    return np.frombuffer(buf, dtype=dtype, count=count, offset=p + 1).copy()


def _load_uncollated(fname, dtype):
    """Read a plain (single-processor) OpenFOAM IOField/IOList."""
    return _parse_block(open(fname, 'rb').read(), dtype)


def _load_collated(fname, dtype, block):
    """Extract one processor block from an OpenFOAM collated (decomposedBlockData)
    file. Blocks are located by their '// processorN <nbytes> (' markers, using
    the byte-length prefixes to skip the (binary) payloads. The master header has
    no such marker and the '// * * *' separator is not always present, so the scan
    starts from the beginning of the file."""
    data = open(fname, 'rb').read()
    pat = re.compile(rb'//\s*processor\d+\s+(\d+)\s*\(')
    pos = 0
    spans = []
    while True:
        mo = pat.search(data, pos)
        if not mo:
            break
        nbytes = int(mo.group(1))
        start = mo.end()
        spans.append((start, nbytes))
        pos = start + nbytes
    start, nbytes = spans[block]
    return _parse_block(data[start:start + nbytes], dtype)


def load(fname, dtype=np.float64, block=None):
    if block is None:
        return _load_uncollated(fname, dtype)
    return _load_collated(fname, dtype, block)


def read_info(fname):
    """Read the nIterations entry from a solve's info dict (the collated info is
    a single global dict, so no per-block handling is required)."""
    info = {}
    with open(fname, errors='replace') as f:
        for raw in f:
            if not raw[0] in {'/', ' '}:
                tok = raw.split()
                if len(tok) == 2:
                    entry = float(tok[1].rstrip(';'))
                    if not tok[0][-4:] in {'ance', 'dual', 'Time'}:
                        entry = int(entry)
                    info[tok[0]] = entry
    return info


def openfoam_residual(A, b, x):

    Ax = A @ x
    Axbar = A @ np.ones(x.shape) * x.mean()
    return np.absolute(b - Ax).sum() / (np.absolute(Ax-Axbar) + np.absolute(b-Axbar)).sum()


def discover(root):
    """Map (field, timestep, solve) -> solve directory for one subdomain root.

    Layout: <root>/<timestep>/<field>-Absol/<solve>/<files>, with an extra
    nesting level (e.g. <timestep>/<region>/<field>-Absol/...) also handled.
    `field` is the path of the dump folder under the timestep with the '-Absol'
    marker dropped: the bare field name ('p_rgh') for a single region, or
    'region/field' ('liquid/p_rgh') when decomposed by region.
    """
    out = {}
    for tdir in glob(f'{root}/*'):
        try:
            timestep = float(tdir.rsplit('/', 1)[-1])
        except ValueError:
            continue
        for depth in (1, 2):
            for fdir in glob(tdir + '/*' * depth):
                if not os.path.isdir(fdir):
                    continue
                comps = fdir.split('/')[-depth:]
                if comps[-1].endswith('-Absol'):
                    comps[-1] = comps[-1][:-len('-Absol')]
                field = '/'.join(comps)
                for sdir in glob(f'{fdir}/*'):
                    try:
                        solve = int(sdir.rsplit('/', 1)[-1])
                    except ValueError:
                        continue
                    out[(field, timestep, solve)] = sdir
    return out


def subdomains(folder):
    """Return [(root, block), ...]. block is an int for a collated processorsN
    directory (one entry per processor block), else None (uncollated processorN
    directories, or a serial / reconstructed tree)."""
    coll = sorted(d for d in glob(f'{folder}/processors*')
                  if re.fullmatch(r'processors\d+', os.path.basename(d)))
    if coll:
        specs = []
        for root in coll:
            n = int(re.fullmatch(r'processors(\d+)', os.path.basename(root)).group(1))
            specs += [(root, b) for b in range(n)]
        return specs
    unc = sorted((d for d in glob(f'{folder}/processor*')
                  if re.fullmatch(r'processor\d+', os.path.basename(d))),
                 key=lambda d: int(os.path.basename(d)[len('processor'):]))
    if unc:
        return [(d, None) for d in unc]
    return [(folder, None)]


def build_orig(specs, region, sizes, offs, nGlobal):
    """Permutation from the decomposed (globalIndex) cell order to the original,
    undecomposed cell order, read from each subdomain's cellProcAddressing. So the
    assembled matrix carries the original mesh numbering and a parallel run matches
    a serial one element-wise. Returns None when the addressing is absent (serial /
    reconstructed tree), i.e. the ordering is already the original one."""
    orig = np.empty(nGlobal, dtype=np.int64)
    for si, (root, block) in enumerate(specs):
        cpa = os.path.join(root, 'constant', *region, 'polyMesh', 'cellProcAddressing')
        if not os.path.isfile(cpa):
            return None
        a = load(cpa, np.uint32, block).astype(np.int64)   # proc-local cell -> original cell
        if a.shape[0] != sizes[si]:
            return None
        orig[offs[si]:offs[si + 1]] = a
    return orig


def desymmetrize(A):
    """Halve the off-diagonal storage of a symmetric matrix: keep the diagonal
    and the doubled strict-upper triangle, so the true matrix is recovered as
    0.5*(A + A.T). This is the form written to the HDF5 matrix datasets."""
    M = A.tocoo()
    dm = M.row == M.col
    um = M.row < M.col
    return coo_array(
        (np.concatenate([M.data[dm], 2.0 * M.data[um]]),
         (np.concatenate([M.row[dm], M.row[um]]),
          np.concatenate([M.col[dm], M.col[um]]))),
        shape=M.shape,
    )


class LinearSystem:
    """One assembled global linear system A x = b. `A` is the symmetric matrix
    (CSR, the true matrix -- not the de-symmetrized storage form); `b`, `sol` and
    `init` are the right-hand side, dumped solution and initial guess; `info` is
    the solver info dict."""

    __slots__ = ('field', 'timestep', 'solve', 'A', 'b', 'sol', 'init', 'info')

    def __init__(self, field, timestep, solve, A, b, sol, init, info):
        self.field, self.timestep, self.solve = field, timestep, solve
        self.A, self.b, self.sol, self.init, self.info = A, b, sol, init, info

    @property
    def niter(self):
        return self.info.get('nIterations', 0)

    def __repr__(self):
        return (f'<LinearSystem {self.field} t={self.timestep} '
                f'solve={self.solve} n={self.A.shape[0]}>')


class LinearSystems:
    """Iterate the linear systems dumped under a case directory, assembling each
    into a monolithic global LinearSystem. Iteration reads from disk but writes
    nothing (run this module's __main__ to persist a whole case to one HDF5 file).

    `folder` is either a case directory (a dump tree) or a single HDF5 file
    written by __main__; the two iterate identically. The HDF5 file holds the
    already-assembled systems with the same carry-forward the tree uses (the
    constant matrix structure stored once, the values only when they change, the
    RHS/solution/initial-guess every solve), and only one system -- at most one
    matrix plus one solve's vectors -- is read into memory at a time.

    By default systems are yielded in the order they were solved: by timestep,
    then by the per-field dump counter. `fields` restricts to the named field(s)
    -- a single name or an iterable, named as `discover` keys them: the bare
    field ('p_rgh') for a single region, or the dump-folder path 'region/field'
    ('liquid/p_rgh') when decomposed by region. `corrector` selects which of a
    timestep's repeated solves to keep: None (default) yields every solve, while
    an integer keeps one solve per (field, timestep), chosen by position among
    that timestep's correctors -- 0 the first, -1 the last."""

    def __init__(self, folder, fields=None, corrector=None):
        self.folder = folder.rstrip('/')
        if fields is None:
            self.fields = None
        elif isinstance(fields, str):
            self.fields = {fields}
        else:
            self.fields = set(fields)
        if not corrector is None and not isinstance(corrector, int):
            raise ValueError("corrector must be None or an int index (e.g. 0 or -1)")
        self.corrector = corrector
        self.h5path = None
        if os.path.isfile(self.folder) and h5py.is_hdf5(self.folder):
            self._init_h5()
        else:
            self.specs = subdomains(self.folder)
            self.maps = {}
            for root, _ in self.specs:
                self.maps.setdefault(root, discover(root))

    def _init_h5(self):
        """Read the per-solve index of an HDF5 dump into memory: the
        (field, timestep, solve) key of each stored system and the ids of its
        carried-forward matrix structure and values. Only this small metadata is
        loaded here; the vectors and matrix arrays are read during iteration."""
        self.h5path = self.folder
        self._h5_index = {}
        with h5py.File(self.h5path, 'r') as f:
            for name, grp in f['systems'].items():
                key = (str(grp.attrs['field']),
                       float(grp.attrs['timestep']),
                       int(grp.attrs['solve']))
                self._h5_index[key] = (name,
                                       int(grp.attrs['struct']),
                                       int(grp.attrs['data']))

    def _ordered_keys(self):
        """Every dumped solve key in solve order (by timestep, then per-field dump
        counter, then field name), restricted to `fields`. This is the full
        carry-forward sequence; corrector selection and completeness are applied
        later in _walk, because a skipped solve still advances the deduplication
        carry-forward and so must stay in the sequence."""
        if self.h5path is not None:
            keys = set(self._h5_index)
        else:
            keys = set().union(*(set(self.maps[r]) for r, _ in self.specs))
        if self.fields is not None:
            keys = {k for k in keys if k[0] in self.fields}
        return sorted(keys, key=lambda k: (k[1], k[2], k[0]))

    def _selected(self, ordered):
        """The keys to emit under `corrector`: all of them for None, else one
        solve per (field, timestep) picked by position among that timestep's
        correctors (`ordered` is ascending in dump counter within each group, so
        index 0 is the first corrector and -1 the last; an out-of-range index
        drops that timestep)."""
        if self.corrector is None:
            return None
        groups = defaultdict(list)
        for key in ordered:
            groups[key[:2]].append(key)          # group by (field, timestep)
        return {grp[self.corrector] for grp in groups.values()
                if -len(grp) <= self.corrector < len(grp)}

    def _walk(self, load_data):
        """Walk every solve in order and yield (key, blocks) for each solve to be
        emitted. The `load_data` argument specifies that `blocks` is empty (if 
        False) or that it carries the emitted solve's carried-forward matrix data 
        (if True)."""
        ordered = self._ordered_keys()
        selected = self._selected(ordered)
        state = defaultdict(dict)
        seen_addr = set()
        for key in ordered:
            field = key[0]
            blocks, complete, addressed = [], True, True
            for si, (root, block) in enumerate(self.specs):
                sd = self.maps[root].get(key)
                if sd is None or not os.path.isfile(f'{sd}/info'):
                    complete = False
                    break
                if load_data:
                    st = state[(si, field)]
                    for q in DEDUP_QUANTS:
                        path = f'{sd}/{q}'
                        if os.path.isfile(path):
                            st[q] = load(path, np.uint32 if q in INT_QUANTS else np.float64, block)
                    blocks.append((dict(st), sd, block))
                    addr = 'lowerAddr' in st
                else:
                    if (si, field) not in seen_addr and os.path.isfile(f'{sd}/lowerAddr'):
                        seen_addr.add((si, field))
                    addr = (si, field) in seen_addr
                if not addr:
                    addressed = False   # addressing not dumped yet (retain the run's first dump)
            if not complete or not addressed:
                continue
            if selected is not None and key not in selected:
                continue   # not the requested corrector: carry-forward advanced, not yielded
            yield key, blocks

    def _walk_h5(self):
        """The HDF5 analogue of _walk: yield (key, (group, struct_id, data_id))
        for each system to emit under `fields`/`corrector`."""
        ordered = self._ordered_keys()
        selected = self._selected(ordered)
        for key in ordered:
            if selected is not None and key not in selected:
                continue
            yield key, self._h5_index[key]

    def keys(self):
        """The (field, timestep, solve) keys that iteration yields, in solve order.
        """
        if self.h5path is not None:
            return [key for key, _ in self._walk_h5()]
        return [key for key, _ in self._walk(load_data=False)]

    def __len__(self):
        if self.h5path is not None:
            return sum(1 for _ in self._walk_h5())
        return sum(1 for _ in self._walk(load_data=False))

    def __iter__(self):
        if self.h5path is not None:
            yield from self._iter_h5()
            return
        orig_cache = {}
        for key, blocks in self._walk(load_data=True):
            yield self._assemble(key, blocks, orig_cache)

    def _iter_h5(self):
        """Yield one assembled system at a time from the HDF5 file. The matrix
        structure and values are deduplicated (carried forward), so consecutive
        systems that reuse them reload nothing: at most one matrix plus one
        solve's b/sol/init vectors are held in memory."""
        with h5py.File(self.h5path, 'r') as f:
            structs, datas, systems = f['matrices/struct'], f['matrices/data'], f['systems']
            si = di = None
            indptr = indices = data = None
            for key, (name, s, d) in self._walk_h5():
                if s != si:
                    g = structs[str(s)]
                    indptr, indices, si = g['indptr'][:], g['indices'][:], s
                if d != di:
                    data, di = datas[str(d)][:], d
                n = indptr.shape[0] - 1
                Ad = csr_array((data, indices, indptr), shape=(n, n))
                A = (0.5 * (Ad + Ad.T)).tocsr()
                grp = systems[name]
                yield LinearSystem(key[0], key[1], key[2], A,
                                   grp['b'][:], grp['sol'][:], grp['init'][:],
                                   json.loads(grp.attrs['info']))

    def _assemble(self, key, blocks, orig_cache):
        """Build the monolithic global system for one (yielded) solve from the
        carried-forward matrix data in `blocks` plus this solve's b/sol/init."""
        field, timestep, solve = key
        info = read_info(f'{blocks[0][1]}/info')

        # --- per-subdomain offsets in the decomposed (globalIndex) order
        sizes = [st['diag'].shape[0] for st, _, _ in blocks]
        offs = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        nGlobal = int(offs[-1])

        # --- map that order back to the original undecomposed numbering via
        #     cellProcAddressing (cached per region; identity when not decomposed)
        region = tuple(os.path.relpath(blocks[0][1], self.specs[0][0]).split('/')[1:-2])
        if region not in orig_cache:
            orig_cache[region] = build_orig(self.specs, region, sizes, offs, nGlobal)
        orig = orig_cache[region]
        if orig is None:
            orig = np.arange(nGlobal, dtype=np.int64)

        # --- assemble the monolithic global symmetric matrix (original ordering),
        #     reading this solve's RHS / solution / initial guess straight from disk
        rows, cols, vals = [], [], []
        bvec = np.zeros(nGlobal)
        solvec = np.zeros(nGlobal)
        initvec = np.zeros(nGlobal)
        for (st, sd, block), lo, hi in zip(blocks, offs[:-1], offs[1:]):
            go = orig[lo:hi]                                             # this block's original cell ids
            rows += [go, orig[st['upperAddr']], orig[st['lowerAddr']]]   # diagonal + internal
            cols += [go, orig[st['lowerAddr']], orig[st['upperAddr']]]
            vals += [st['diag'], st['lower'], st['lower']]
            ir = st.get('interfaceRow')                                  # interface couplings
            if ir is not None and ir.size:
                rows.append(orig[ir]); cols.append(orig[st['interfaceCol']]); vals.append(st['interfaceLower'])
            bvec[go] = load(f'{sd}/b', np.float64, block)
            solvec[go] = load(f'{sd}/sol', np.float64, block)
            try:
                init_b = load(f'{sd}/init', np.float64, block)
            except (ValueError, OSError):
                init_b = None
            if init_b is not None and init_b.size:
                initvec[go] = init_b

        A = coo_array(
            (np.concatenate(vals).astype(np.float64),
             (np.concatenate(rows).astype(np.int64),
              np.concatenate(cols).astype(np.int64))),
            shape=(nGlobal, nGlobal),
        ).tocsr()   # tocsr() sums any coincident entries; this is the true symmetric matrix
        return LinearSystem(field, timestep, solve, A, bvec, solvec, initvec, info)


if __name__ == '__main__':

    # Repackage a whole dump tree into one HDF5 file that LinearSystems iterates
    # exactly as it does the tree. Usage: Absol.py <case_dir> [out.h5].
    folder = sys.argv[1].rstrip('/')
    out = sys.argv[2] if len(sys.argv) > 2 else folder + '.h5'

    systems = LinearSystems(folder)
    total = len(systems)

    mrr = 0.0
    ml2 = 0.0
    struct_state = {}    # field -> (id, indptr, indices) most recently stored
    data_state = {}      # field -> (id, values) most recently stored
    si = di = -1         # next structure / values id to allocate

    with h5py.File(out, 'w') as h5:
        gsys = h5.create_group('systems')
        gstruct = h5.create_group('matrices/struct')
        gdata = h5.create_group('matrices/data')

        for n, system in enumerate(systems):

            # De-symmetrized (diag + doubled strict upper) CSR storage. The
            # structure (indptr/indices) is constant for a field over the whole
            # run and its values change only once per timestep, so each is
            # written only when it differs from that field's previous solve.

            field = system.field
            Ad = desymmetrize(system.A).tocsr()

            ss = struct_state.get(field)
            if ss is None or not (np.array_equal(Ad.indptr, ss[1])
                                  and np.array_equal(Ad.indices, ss[2])):
                si += 1
                g = gstruct.create_group(str(si))
                g.create_dataset('indptr', data=Ad.indptr)
                g.create_dataset('indices', data=Ad.indices)
                struct_state[field] = ss = (si, Ad.indptr, Ad.indices)

            ds = data_state.get(field)
            if ds is None or not np.array_equal(Ad.data, ds[1]):
                di += 1
                gdata.create_dataset(str(di), data=Ad.data)
                data_state[field] = ds = (di, Ad.data)

            gs = gsys.create_group(str(n))
            gs.create_dataset('b', data=system.b)
            gs.create_dataset('sol', data=system.sol)
            gs.create_dataset('init', data=system.init)
            gs.attrs['field'] = system.field
            gs.attrs['timestep'] = system.timestep
            gs.attrs['solve'] = system.solve
            gs.attrs['struct'] = ss[0]
            gs.attrs['data'] = ds[0]
            gs.attrs['info'] = json.dumps(system.info)

            print('processed', n + 1, '/', total, 'systems', end='\t')
            effTol = max(system.info['relativeTolerance'] * openfoam_residual(system.A, system.b, system.init),
                         system.info['tolerance'])
            mrr = max(mrr, openfoam_residual(system.A, system.b, system.sol) / effTol)
            print('max residual:reported:', round(mrr, 2), end='\t')
            ml2 = max(ml2, np.linalg.norm(system.A @ system.sol - system.b) / np.linalg.norm(system.b))
            print('max l2 error:', ml2, end='\r')

    print()
