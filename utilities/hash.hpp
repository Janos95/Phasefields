//
// Created by janos on 18.04.20.
//

#pragma once

//inline size_t hash_int(uint64_t obj) noexcept {
//#if ROBIN_HOOD(HAS_UMUL128)
//    // 167079903232 masksum, 120428523 ops best: 0xde5fb9d2630458e9
//    static constexpr uint64_t k = UINT64_C(0xde5fb9d2630458e9);
//    uint64_t h;
//    uint64_t l = detail::umul128(obj, k, &h);
//    return h + l;
//#elif ROBIN_HOOD(BITNESS) == 32
//    uint64_t const r = obj * UINT64_C(0xca4bcaa75ec3f625);
//    auto h = static_cast<uint32_t>(r >> 32U);
//    auto l = static_cast<uint32_t>(r);
//    return h + l;
//#else
//    // murmurhash 3 finalizer
//    uint64_t h = obj;
//    h ^= h >> 33;
//    h *= 0xff51afd7ed558ccd;
//    h ^= h >> 33;
//    h *= 0xc4ceb9fe1a85ec53;
//    h ^= h >> 33;
//    return static_cast<size_t>(h);
//#endif
//}


