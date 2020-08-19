//
// taken from llvm STLExtras.h
//

#pragma once

#include <type_traits>
#include <cstdint>

/// An efficient, type-erasing, non-owning reference to a callable. This is
/// intended for use as the type of a function parameter that is not used
/// after the function in question returns.
///
/// This class does not own the callable, so it is not in general safe to store
/// a function_ref.
template<typename Fn>
class FunctionRef;

template<typename Ret, typename ...Params>
class FunctionRef<Ret(Params...)> {
    Ret (* callback)(intptr_t callable, Params ...params) = nullptr;

    intptr_t callable;

    template<typename Callable>
    static Ret callback_fn(intptr_t callable, Params ...params) {
        return (*reinterpret_cast<Callable*>(callable))(
                (Params&&) (params)...);
    }

public:
    FunctionRef() = default;

    FunctionRef(std::nullptr_t) {}

    template<typename Callable>
    FunctionRef(Callable&& callable,
                typename std::enable_if<
                        !std::is_same<typename std::remove_reference<Callable>::type,
                                FunctionRef>::value>::type* = nullptr)
            : callback(callback_fn<typename std::remove_reference<Callable>::type>),
              callable(reinterpret_cast<intptr_t>(&callable)) {}

    Ret operator()(Params ...params) const {
        return callback(callable, (Params&&) (params)...);
    }

    explicit operator bool() const { return callback; }
};

