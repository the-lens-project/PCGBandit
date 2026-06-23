#!/bin/bash

docker pull opencfd/openfoam-default:latest

docker run -it --rm \
  -v "$(pwd)":/home/openfoam \
  opencfd/openfoam-default:latest /bin/bash -lc '
set -e

cat > /tmp/fakeroot.c << "EOF"
#include <sys/types.h>
uid_t geteuid(void) { return 1000; }
EOF

gcc -shared -o /tmp/fakeroot.so /tmp/fakeroot.c
export LD_PRELOAD=/tmp/fakeroot.so

cd /home/openfoam/src/PCGBandit
wmake libso
cd ../ICTC
wmake libso
cd ../FGAMG
wmake libso
cd /home/openfoam

exec /bin/bash
'
