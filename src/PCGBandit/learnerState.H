/*---------------------------------------------------------------------------*\
Class
    Foam::learnerState

Description
    In-memory, move-only learner state shared between solver instances.
    Unlike dictionary, keys are literal and values retain their exact C++
    types; lookups return references without parsing or serialization.
    References survive insertion, in-place updates, and ownership transfer.
    Replacing a state destroys its previous entries; assigning a field may
    reallocate its element buffer. Reusing a key with another type is an error.
\*---------------------------------------------------------------------------*/

#ifndef learnerState_H
#define learnerState_H

#include "HashPtrTable.H"
#include "autoPtr.H"
#include "error.H"
#include "word.H"
#include "typeInfo.H"
#include <memory>
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace Foam
{

class learnerState
{
    struct entry
    {
        //- Destroy the typed payload through the common entry base.
        virtual ~entry() = default;
    };

    template<class T>
    struct storedValue : entry
    {
        T value;

        //- Construct T directly from arguments, without a token stream.
        template<class... Args>
        explicit storedValue(Args&&... args)
        :
            value(std::forward<Args>(args)...)
        {}
    };

    // The table owns each polymorphic entry; rehashing preserves references.
    HashPtrTable<entry> entries_;

    //- Insert or assign the same type in place, skipping self-assignment.
    template<class T, class Value>
    T& setValue(const word& keyword, Value&& value)
    {
        if (!found(keyword))
        {
            return emplace<T>(keyword, std::forward<Value>(value));
        }
        T& existing = get<T>(keyword);
        // In particular, avoid OpenFOAM Field self-assignment.
        if (std::addressof(existing) != std::addressof(static_cast<const T&>(value)))
        {
            existing = std::forward<Value>(value);
        }
        return existing;
    }

public:

    //- Same as dictionary: default-construct an empty container.
    learnerState() = default;

    //- Different from dictionary: copying is disabled; payloads need not copy.
    learnerState(const learnerState&) = delete;

    //- Different from dictionary: copy assignment is disabled.
    learnerState& operator=(const learnerState&) = delete;

    //- Same as dictionary: move-construct by transferring entry ownership.
    learnerState(learnerState&&) = default;

    //- Different from dictionary: move assignment transfers entries, not copies.
    learnerState& operator=(learnerState&&) = default;

    //- Different from dictionary::found: literal local keys only, no patterns.
    bool found(const word& keyword) const
    {
        return entries_.found(keyword);
    }

    //- Different from dictionary::get: return const T& of the exact stored type.
    template<class T>
    const T& get(const word& keyword) const
    {
        const entry* value = entries_.get(keyword);
        if (!value)
        {
            FatalErrorInFunction
                << "Missing learner state entry '" << keyword << "'"
                << abort(FatalError);
        }
        const auto* stored = isA<storedValue<T>>(*value);
        if (!stored)
        {
            FatalErrorInFunction
                << "Learner state entry '" << keyword
                << "' does not hold the requested type " << typeid(T).name()
                << abort(FatalError);
        }
        return stored->value;
    }

    //- Different from dictionary::get: return mutable T&, not a parsed value.
    template<class T>
    T& get(const word& keyword)
    {
        return const_cast<T&>(static_cast<const learnerState&>(*this).get<T>(keyword));
    }

    //- No dictionary counterpart: construct T and return T&; duplicate keys are fatal.
    template<class T, class... Args>
    T& emplace(const word& keyword, Args&&... args)
    {
        if (found(keyword))
        {
            FatalErrorInFunction
                << "Learner state entry '" << keyword << "' already exists"
                << abort(FatalError);
        }
        entries_.insert
        (
            keyword,
            autoPtr<entry>::NewFrom<storedValue<T>>(std::forward<Args>(args)...)
        );
        return get<T>(keyword);
    }

    //- Different from dictionary::set: insert/assign same-typed T and return T&.
    template<class T>
    T& set(const word& keyword, const T& value)
    {
        return setValue<T>(keyword, value);
    }

    //- Different from dictionary::set: copy/move decayed T in place; return T&.
    template<class T>
    auto& set(const word& keyword, T&& value)
    {
        return setValue<std::decay_t<T>>
        (
            keyword, std::forward<T>(value)
        );
    }

    //- Different from dictionary::getOrAdd: return T&; copy default if absent.
    template<class T>
    T& getOrAdd(const word& keyword, const T& defaultValue)
    {
        return found(keyword) ? get<T>(keyword) : emplace<T>(keyword, defaultValue);
    }

    //- Different from dictionary::getOrAdd: return T&; copy/move default if absent.
    template<class T>
    auto& getOrAdd(const word& keyword, T&& defaultValue)
    {
        using valueType = std::decay_t<T>;
        return found(keyword)
            ? get<valueType>(keyword)
            : emplace<valueType>(keyword, std::forward<T>(defaultValue));
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
