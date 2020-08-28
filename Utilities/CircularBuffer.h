//
// Created by janos on 8/18/20.
//

#pragma once

#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield::Containers {

namespace Cr = Corrade;

namespace {

std::uint32_t nextPow2(std::uint32_t v) {
    v--;
    v |= v >> 1u;
    v |= v >> 2u;
    v |= v >> 4u;
    v |= v >> 8u;
    v |= v >> 16u;
    v++;
    return v;
}

}

template<class T>
class CircularBuffer {
public:

    enum class RelocationOption : std::uint8_t {
        LeftBound = 1,
        RightBound = 2,
        Center = 3
    };

    explicit CircularBuffer(std::size_t size) : m_data{Cr::Containers::DefaultInit, nextPow2(size)} {
        m_mask = m_data.size() - 1;
    }

    explicit CircularBuffer() = default;

    template<class... Args>
    void emplaceBack(Args&& ... args) {
        if(m_size + 1 > m_data.size()){
            grow(RelocationOption::Center);
        }
        m_data[m_end] = T{static_cast<Args&&>(args)...};
        m_end = m_mask & (m_end + 1);
        ++m_size;
    }

    template<class... Args>
    void emplaceFront(Args&& ... args) {
        if(m_size + 1 > m_data.size()){
            grow(RelocationOption::Center);
        }
        m_begin = (m_begin - 1) & m_mask;
        m_data[m_begin] = T{static_cast<Args&&>(args)...};
        ++m_size;
    }

    T popBack() {
        --m_size;
        m_end = m_mask & (m_end - 1);
        return std::move(m_data[m_end]);
    }

    T popFront() {
        --m_size;
        auto oldBegin = m_begin;
        m_begin = m_mask & (m_begin + 1);
        return std::move(m_data[oldBegin]);
    }

    void shrinkToFit() {
        auto size = nextPow2(m_size);
        if(size != m_data.size()) {
            Cr::Containers::Array<T> data{Corrade::Containers::DefaultInit, size};

            for(std::uint32_t i = 0; i < m_size; ++i)
                data[i] = std::move(m_data[(m_begin + i) & m_mask]);

            m_data = data;
            m_begin = 0;
            m_end = m_size;
        }
    }

    T& operator[](std::size_t idx) {
        return m_data[(m_begin + idx) & m_mask];
    }

    const T& operator[](std::size_t idx) const {
        return m_data[(m_begin + idx) & m_mask];
    }

    void grow(RelocationOption option = RelocationOption::Center){
        Cr::Containers::Array<T> data(Cr::Containers::DefaultInit, nextPow2(m_size + 1));

        std::uint32_t offset;
        switch(option){
            case RelocationOption::LeftBound : offset = 0; break;
            case RelocationOption::Center : offset = data.size()/2 - m_size/2; break;
            case RelocationOption::RightBound : offset = data.size() - m_size; break;
        }
        for(std::uint32_t i = 0; i < m_size; ++i)
            data[i + offset] = m_data[(m_begin + i) & m_mask];

        m_data = std::move(data);
        m_begin = offset;
        m_end = m_size + offset;
        m_mask = m_data.size() - 1;
    }

    [[nodiscard]] std::uint32_t size() const {
        return m_size;
    }

    [[nodiscard]] bool empty() const {
        return m_size == 0;
    }

private:
    std::uint32_t m_begin = 0, m_end = 0, m_size = 0;
    Cr::Containers::Array<T> m_data;

    std::uint32_t m_mask = 0;
};


}
