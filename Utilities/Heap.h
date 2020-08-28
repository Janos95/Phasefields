//
// Created by janos on 8/18/20.
//

#pragma once

#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield::Containers {

namespace Implementation {

struct Comparator {
    template<class T>
    bool operator()(T const& x, T const& y) {
        return x < y;
    }
};

}

template<class T, class Comp = Implementation::Comparator>
class Heap {
public:

    explicit Heap(Comp comp = {}) {}

    T const& findMin() const {
        CORRADE_ASSERT(!empty(), "Heap : find_min failed because Heap is empty", {});
        return _data[0];
    }

    T extractMin() {
        T result = findMin();
        remove(0);
        return result;
    }

    template<class... Args>
    int emplace(Args&& ... args) {
        arrayAppend(_data, Corrade::Containers::InPlaceInit, static_cast<Args&&>(args)...);
        siftUp(_data.size() - 1);
        return _data.size() - 1;
    }

    void clear() {
        Corrade::Containers::arrayResize(_data, 0);
    }

    bool empty() const {
        return _data.empty();
    }

protected:

    void remove(int index) {
        CORRADE_ASSERT(index < _data.size(), "Index error in heap",);
        using std::swap;
        swap(_data[index], _data[_data.size() - 1]);
        arrayResize(_data, _data.size() - 1);
        siftUp(index);
        siftDown(index);
    }

    void descreaseKey(int index) {
        CORRADE_ASSERT(index < _data.size(), "Index error in heap",);
        siftUp(index);
    }

    T& operator[](int index) {
        return _data[index];
    }

    T const& operator[](int index) const {
        return _data[index];
    }

private:

    static std::size_t parent(int index) {
        return (index - 1)/2;
    }

    static std::size_t left(std::size_t index) {
        return (2*index) + 1;
    }

    static std::size_t right(std::size_t index) {
        return (2*index) + 2;
    }

    void siftUp(int index) {
        while((index > 0) and (_data[index] < _data[parent(index)])) {
            using std::swap; /* two phase lookup */
            swap(_data[index], _data[parent(index)]);
            index = parent(index);
        }
    }

    void siftDown(int index) {
        std::size_t smallest = index;
        while(true) {
            if(left(index) <_data.size() && _comp(_data[left(index)], _data[smallest]) {
                smallest = left(index);
            }
            if(right(index) < _data.size() && _data[right(index)] < _data[smallest]) {
                smallest = right(index);
            }
            if(index == smallest) return;
            swap(_data[smallest], _data[index]);
            index = smallest;
        }
    }

    Corrade::Containers::Array<T> _data;       // holds the objects in heap order
    Comp _comp;
};

}