#pragma once

#include <type_traits>

template<typename Fn>
class FunctionRef;

template<typename Ret, typename ...Params>
class FunctionRef<Ret(Params...)> {
    Ret (* callback)(void* callable, Params ...params) = nullptr;
    void* callable;

public:
    FunctionRef() = default;

    FunctionRef(std::nullptr_t) {}

    template<typename Callable>
    FunctionRef(Callable&& callable) : callable(&callable) {
        callback = +[](void* c, Params... params) { return (*static_cast<std::remove_reference_t<Callable>*>(c))(params...); };
    }

    Ret operator()(Params ...params) const {
        return callback(callable, (Params&&) (params)...);
    }

    explicit operator bool() const { return callback; }
};

