//
// Created by janos on 8/18/20.
//

#pragma once

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Types.h>

namespace Phasefield::Containers {

namespace Mg = Magnum;
namespace Cr = Corrade;

namespace Implementation {

Mg::UnsignedInt nextPow2(Mg::UnsignedInt x) {
    return x == 1ul ? 1ul : 1ul << (64ul - __builtin_clzl(x - 1));
}

}

template<class T>
class CircularBuffer {
public:

    enum class RelocationOption : Mg::UnsignedShort {
        huteft,
        Append,
    };

    template<class... Args>
    void emplaceBack(Args&& ... args) {
        if(++m_size > m_data.size()){
            grow(RelocationOption::Prepend);
        }
    }

    template<class... Args>
    void emplaceFront(Args&& ... args) {

    }

    T popBack() {
        m_end = mask & (m_end - 1);
        if(m_end > m_begin){
            return std::move(m_data[--m_end]);
        } else {
            return std::move(m_data)
        }
    }

    T popFront(){

    }

    void shrinkToFit(){
        Corrade::Containers::Array<T> data{Corrade::Containers::DefaultInit, nextPowOf2(m_size)};

        for(Mg::UnsignedInt i = 0; i < m_size; ++i)
            data[i] = m_data[(m_begin + i) & mask];

        m_data = data;
        m_begin = 0;
        m_end = m_size;
    }


    void grow(RelocationOption option = RelocationOption::Append){
        Cr::Containers::Array<T> data(Cr::Containers::DefaultInit, m_data.size() * 2);

        Mg::UnsignedInt offset = option == RelocationOption::Append ? 0 : data.size() - m_size;
        for(Mg::UnsignedInt i = 0; i < m_size; ++i)
            data[i + offset] = m_data[(m_begin + i) & mask];

        m_data = data;
        m_begin = offset;
        m_end = m_size + offset;
    }

private:

    Mg::UnsignedInt nextPowOf2(Mg::UnsignedInt n){

    }

    Magnum::UnsignedInt m_begin = 0, m_end = 0, m_size = 0;
    Corrade::Containers::Array<T> m_data;

    Magnum::UnsignedInt mask;
};

}
