//
// Created by janos on 8/18/20.
//

#pragma once

#include <cstdlib>
#include <cstring>
#include <type_traits>

#ifdef WITH_STD_HEAP
#include <algorithm>
#endif

#include "MinMaxAndDAryHeap.h"

struct MinHeap {
    template<class T>
    bool operator()(T const& x, T const& y) const {
        return x > y;
    }
};

struct MaxHeap {
    template<class T>
    bool operator()(T const& x, T const& y) const {
        return x < y;
    }
};

template<class Comp>
struct SwapComp {
    Comp comp;
    template<class T>
    bool operator()(T const& x, T const& y) const {
        return !comp(x, y);
    }
};

enum class HeapAlgorithm : uint8_t {
    MinMax = 1,
    TwoAry = 2,
    FourAry = 3
#ifdef WITH_STD_HEAP
    ,StdHeap = 4
#endif
};

template<class T, class Comp = MinHeap, HeapAlgorithm algo = HeapAlgorithm::FourAry>
class Heap {
public:

    constexpr static size_t Arity = algo == HeapAlgorithm::TwoAry ? 2 : algo == HeapAlgorithm::FourAry ? 4 : 0;

    explicit Heap(Comp comp = {}) : m_comp{comp} {}

    explicit Heap(size_t n, Comp comp = {}) : m_comp{comp} { reserve(n); }

    template<class It>
    explicit Heap(It first, It end, Comp comp = {}) : m_comp{comp} {
        m_size = end - first;
        m_capacity = m_size;
        m_data = (T*)malloc(sizeof(T)*m_size);

        T* p = m_data;
        for(auto it = first; it != end; ++it, ++p)
            ::new(static_cast<void*>(p)) T{*it};

        if constexpr (algo == HeapAlgorithm::MinMax)
            make_minmax_heap(m_data, m_data + m_size, SwapComp<Comp>{m_comp});
#ifdef WITH_STD_HEAP
        else if constexpr (algo == HeapAlgorithm::StdHeap)
            std::make_heap(m_data, m_data + m_size, m_comp);
#endif
        else
            make_dary_heap<Arity>(m_data, m_data + m_size, m_comp);
    }

    template<class It>
    void insert(It first, It end) {
        bool isEmpty = m_size == 0;
        m_size += end - first;
        reserve(m_size);
        if(isEmpty) { /* if the heap is empty we call the linear subroutine */
            T* p = m_data;
            for(; first != end; ++first, ++p)
                ::new(static_cast<void*>(p)) T{*first};

            if constexpr (algo == HeapAlgorithm::MinMax)
                make_minmax_heap(m_data, m_data + m_size, SwapComp<Comp>{m_comp});
#ifdef WITH_STD_HEAP
            else if constexpr (algo == HeapAlgorithm::StdHeap)
                std::make_heap(m_data, m_data + m_size, m_comp);
#endif
            else
                make_dary_heap<Arity>(m_data, m_data + m_size, m_comp);
        } else {
            for(; first != end ; ++first)
                emplace(*first);
        }
    }

    T const& findMin() const {
        return m_data[0];
    }

    T extractMin() {
        if constexpr (Arity > 0)
            pop_dary_heap<Arity>(m_data, m_data + m_size--, m_comp);
#ifdef WITH_STD_HEAP
        else if constexpr (algo == HeapAlgorithm::StdHeap) {
            std::pop_heap(m_data, m_data + m_size--, m_comp);
        }
#endif
        else pop_minmax_heap_min(m_data, m_data + m_size--, SwapComp<Comp>{m_comp});

        return static_cast<T&&>(m_data[m_size]);
    }

    T const& findMax() const requires (algo == HeapAlgorithm::MinMax) {
        if (m_size <= 2ul)
            return m_data[m_size - 1];
        size_t index = 1ul + static_cast<size_t>(m_comp(m_data[1], m_data[2]));
        return m_data[index];
    }

    T extractMax() requires (algo == HeapAlgorithm::MinMax) {
        pop_minmax_heap_max(m_data, m_data + m_size--, SwapComp<Comp>{m_comp});
        return static_cast<T&&>(m_data[m_size]);
    }

    template<class... Args>
    void emplace(Args&& ... args) {
        if(m_capacity == m_size)
            reserve(computeEmplaceCapacity());

        ::new(m_data + m_size) T{args...};
        ++m_size;
        if constexpr (Arity > 0)
            push_dary_heap<Arity>(m_data, m_data + m_size, m_comp);
#ifdef WITH_STD_HEAP
        else if constexpr (algo == HeapAlgorithm::StdHeap) {
            std::push_heap(m_data, m_data + m_size, m_comp);
        }
#endif
        else push_minmax_heap(m_data, m_data + m_size, SwapComp<Comp>{m_comp});
    }

    T dummy(T x){
        return 2*x;
    }

    void clear() {
        destruct();
        m_size = 0;
    }

    void shrinkToFit() {
        reserve(m_size);
    }

    void reserve(size_t targetCapacity) {
        if(targetCapacity > m_capacity) {
            m_capacity = targetCapacity;
            if constexpr (std::is_trivially_copyable_v<T>) {
                m_data = (T*)std::realloc(m_data, m_capacity*sizeof(T)); /* if m_data is nullptr, realloc behaves like malloc */
            } else { /* we assume nothing throws */
                T* data = std::malloc(m_capacity*sizeof(T));
                T* to = data, from = m_data;
                for(; from != m_data + m_size; ++to, ++from) {
                    ::new(static_cast<void*>(to)) T{(T && )*from};
                    from->~T();
                }
                free(m_data);
                m_data = data;
            }
        }
    }

    [[nodiscard]] bool empty() const { return m_size == 0; }

    T& operator[](int index) { return m_data[index]; }

    T const& operator[](int index) const { return m_data[index]; }

    [[nodiscard]] std::size_t size() const { return m_size; }

    [[nodiscard]] T const* begin() const { return m_data; }

    [[nodiscard]] T const* end() const { return m_data + m_size; }

    [[nodiscard]] T const* data() const { return m_data; }

    ~Heap() {
        destruct();
        free(m_data);
    }

private:

    /* based on folly's vector class */
    size_t computeEmplaceCapacity() const {
        static constexpr size_t mallocMinInPlaceExpandable = 4096;

        if (m_capacity == 0) {
            return heap_detail::max(64/sizeof(T), 1ul);
        }
        if (m_capacity < mallocMinInPlaceExpandable / sizeof(T)) {
            return m_capacity * 2;
        }
        if (m_capacity > mallocMinInPlaceExpandable * 32 / sizeof(T)) {
            return m_capacity * 2;
        }
        return (m_capacity * 3 + 1) / 2;
    }

    void destruct() {
        if constexpr (!std::is_trivially_destructible_v<T>) {
            for(size_t i = 0; i < m_size; ++i) {
                m_data[i].~T();
            }
        }
    }

    T* m_data = nullptr;
    size_t m_size = 0;
    size_t m_capacity = 0;
    Comp m_comp;
};
