//
// Created by janos on 12.04.20.
//

#pragma once

#include <Corrade/Containers/StaticArray.h>

namespace Cr = Corrade;

template<int N, class T>
class SmallArray {
public:
    using value_type = T;
    using reference = T&;
    using const_reference = T const&;
    using iterator = T*;
    using const_iterator = T const*;
    using difference_type = std::ptrdiff_t;
    using size_type = int;


    explicit SmallArray(Cr::Containers::DefaultInitT) : m_data{Cr::Containers::DefaultInit} {}

    explicit SmallArray(Cr::Containers::ValueInitT) : m_data{Cr::Containers::ValueInit} {}

    explicit SmallArray(Cr::Containers::NoInitT) : m_data(Cr::Containers::NoInit) {}

    iterator begin() { return m_data.data(); }

    const_iterator begin() const { return m_data.data(); }

    iterator end() { return m_data.data() + m_size; }

    const_iterator end() const { return m_data.data() + m_size; }

    reference operator[](int i) {
        return m_data[i];
    }

    const_reference operator[](int i) const {
        return m_data[i];
    }

    template<class... Args>
    void emplace_back(Args&& ... args) {
        m_data[m_size++] = T((Args&&) args...);
    }

    void push_back(T&& x) {
        m_data[m_size++] = (T&&) x;
    }

    [[nodiscard]] size_type size() const { return m_size; }

private:

    Corrade::Containers::StaticArray<N, T> m_data;
    std::size_t m_size = 0;
};