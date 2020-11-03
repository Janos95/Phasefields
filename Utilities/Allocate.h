//
// Created by janos on 10/10/20.
//

#pragma once

#include <new>

namespace Phasefield {

inline void* allocate_buffer(size_t Size, size_t Alignment) {
    return ::operator new(Size
#ifdef __cpp_aligned_new
            ,
                          std::align_val_t(Alignment)
#endif
    );
}


inline void deallocate_buffer(void* Ptr, size_t Size, size_t Alignment) {
    ::operator delete(Ptr
#ifdef __cpp_sized_deallocation
            ,
                      Size
#endif
#ifdef __cpp_aligned_new
            ,
                      std::align_val_t(Alignment)
#endif
    );
}



}

