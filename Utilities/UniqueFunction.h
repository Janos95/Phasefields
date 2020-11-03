#pragma once

#include "Allocate.h"
#include "Utility.h"

namespace Phasefield {

template<typename FunctionT>
class UniqueFunction;

template<typename Ret, typename ...Params>
class UniqueFunction<Ret(Params...)> {
public:

    UniqueFunction() = default;

    UniqueFunction(std::nullptr_t) {}

    template<class Callable>
    UniqueFunction(Callable&& f) {

        erased = allocate_buffer(sizeof(Callable), alignof(Callable));
        ::new(erased) Callable((Callable&&) f);

        destroy = +[](void* e) {
            static_cast<Callable*>(e)->~Callable();
            deallocate_buffer(e, sizeof(Callable), alignof(Callable));
        };

        call = +[](void* e, Params... params) -> Ret {
            return (*static_cast<Callable*>(e))(params...);
        };
    }

    UniqueFunction& operator=(UniqueFunction&& other) noexcept {
        other.swap(*this);
        return *this;
    }

    UniqueFunction(UniqueFunction&& other) noexcept {
        other.swap(*this);
    }

    void swap(UniqueFunction& other) {
        pf_swap(erased, other.erased);
        pf_swap(destroy, other.destroy);
        pf_swap(call, other.call);
    }

    Ret operator()(Params ...params) const {
        return call(erased, (Params&&) (params)...);
    }

    explicit operator bool() const { return erased; }

    ~UniqueFunction() {
        if(destroy) destroy(erased);
    }

private:

    void* erased = nullptr;
    void (* destroy)(void*) = nullptr;
    Ret (* call)(void*, Params...) = nullptr;
};

}
