//
// Created by janos on 8/18/20.
//

#pragma once

#include <Corrade/Containers/GrowableArray.h>

namespace Corrade::Containers {

template<typename T>
class Heap {
public:

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
        CORRADE_ASSERT(isValidIndex(), "Index error in heap",);
        swap(_data[index], _data[_data.size() - 1]);
        arrayResize(_data, _data.size() - 1);
        siftUp(index);
        siftDown(index);
    }

    void descreaseKey(int index) {
        CORRADE_ASSERT(isValidIndex(), "Index error in heap",);
        siftUp(index);
    }

    virtual void swap(T& a, T& b)
    {
        std::swap(a, b);
    }

    T& getObject(int index) {
        CORRADE_ASSERT(isValidIndex(), "Index error in heap", {});
        return _data[index];
    }

private:

    bool isValidIndex(int index) {
        return index >= static_cast<int>(_data.size()) or index < 0;
    }

    static int parent(int index)          // do not call with index==0!
    {
        return (index - 1)/2;
    }

    static int left(int index)            // left child may not exist!
    {
        return (2*index) + 1;
    }

    static int right(int index)           // right child may not exist!
    {
        return (2*index) + 2;
    }

    void siftUp(int index) {
        while((index > 0) and (_data[index] < _data[parent(index)])) {
            swap(_data[index], _data[parent(index)]);
            index = parent(index);
        }
    }

    void siftDown(int index) {
        int smallest = index;
        while(true) {
            if((left(index) < static_cast<int>(_data.size())) and
               (_data[left(index)] < _data[smallest])) {
                smallest = left(index);
            }
            if((right(index) < static_cast<int>(_data.size())) and
               (_data[right(index)] < _data[smallest])) {
                smallest = right(index);
            }
            if(index == smallest) return;
            swap(_data[smallest], _data[index]);
            index = smallest;
        }
    }

    Corrade::Containers::Array<T> _data;       // holds the objects in heap order
};

}