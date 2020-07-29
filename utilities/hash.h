//
// Created by janos on 18.04.20.
//

#pragma once

#include <cstdint>
#include <cstddef>

#define ROBIN_HOOD(x) ROBIN_HOOD_PRIVATE_DEFINITION_##x()

namespace detail {

// umul
#if defined(__SIZEOF_INT128__)
#    define ROBIN_HOOD_PRIVATE_DEFINITION_HAS_UMUL128() 1
#    if defined(__GNUC__) || defined(__clang__)
#        pragma GCC diagnostic push
#        pragma GCC diagnostic ignored "-Wpedantic"
using uint128_t = unsigned __int128;
#        pragma GCC diagnostic pop
#    endif

inline uint64_t umul128(uint64_t a, uint64_t b, uint64_t* high) noexcept {
    auto result = static_cast<uint128_t>(a)*static_cast<uint128_t>(b);
    *high = static_cast<uint64_t>(result >> 64U);
    return static_cast<uint64_t>(result);
}

#elif (defined(_MSC_VER) && ROBIN_HOOD(BITNESS) == 64)
#    define ROBIN_HOOD_PRIVATE_DEFINITION_HAS_UMUL128() 1
#    include <intrin.h> // for __umulh
#    pragma intrinsic(__umulh)
#    ifndef _M_ARM64
#        pragma intrinsic(_umul128)
#    endif
inline uint64_t umul128(uint64_t a, uint64_t b, uint64_t* high) noexcept {
#    ifdef _M_ARM64
    *high = __umulh(a, b);
    return ((uint64_t)(a)) * (b);
#    else
    return _umul128(a, b, high);
#    endif
}
#else
#    define ROBIN_HOOD_PRIVATE_DEFINITION_HAS_UMUL128() 0
#endif
}

inline size_t hash_int(uint64_t obj) noexcept {
#if ROBIN_HOOD(HAS_UMUL128)
    // 167079903232 masksum, 120428523 ops best: 0xde5fb9d2630458e9
    static constexpr uint64_t k = UINT64_C(0xde5fb9d2630458e9);
    uint64_t h;
    uint64_t l = detail::umul128(obj, k, &h);
    return h + l;
#elif ROBIN_HOOD(BITNESS) == 32
    uint64_t const r = obj * UINT64_C(0xca4bcaa75ec3f625);
    auto h = static_cast<uint32_t>(r >> 32U);
    auto l = static_cast<uint32_t>(r);
    return h + l;
#else
    // murmurhash 3 finalizer
    uint64_t h = obj;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccd;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53;
    h ^= h >> 33;
    return static_cast<size_t>(h);
#endif
}

