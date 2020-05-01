//
// Created by janos on 29.04.20.
//

#pragma once

#include <type_traits> //for is_trivial
#include <utility> //for move
#include <cstring> //for memcpy
#include <cstdlib> //for malloc

template<class T>
struct SharedRessource {

    static_assert(std::is_trivial_v<T>, "Shared Ressource : T needs to be trivial");
    struct Block { T x; size_t refCount; };
    Block * data = nullptr;

    SharedRessource(nullptr_t) : data(nullptr) {}

    explicit SharedRessource(T const& x)
    {
        data = (Block*)std::malloc(sizeof(Block));
        std::memcpy(&data->x, &x, sizeof(T));
        data->refCount = 1;
    }

    SharedRessource& operator=(SharedRessource other) noexcept {
        other.swap(*this);
        return *this;
    }

    SharedRessource(SharedRessource const& other) noexcept : data(other.data) {
        ++(data->refCount);
    }

    SharedRessource(SharedRessource&& other) noexcept : data(other.data) {
        other.data = nullptr;
    }

    void swap(SharedRessource& other){
        auto temp = other.data;
        other.data = data;
        data = temp;
    }

    T& operator*(){ return data->x; }
    T const& operator*() const { return data->x; }

    explicit operator bool() const { return data != nullptr; }

    T* get() { return &data->x; }
    T const* get() const { return &data->x; }

    int refCount() { return data->refCount; }

    ~SharedRessource(){
        if(!data) return;
        CORRADE_INTERNAL_ASSERT(data->refCount);
        if(!(--(data->refCount)))
            delete data;
    }
};
