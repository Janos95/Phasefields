//
// Created by janos on 9/9/20.
//

#ifndef PHASEFIELD_SERIALIZE_H
#define PHASEFIELD_SERIALIZE_H

#include "Types.h"
#include <Corrade/Containers/GrowableArray.h>

template<class T>
void serializeTrivial(Array<char>& data, T const& x) {
    arrayAppend(data, {&x, sizeof(std::remove_reference_t<T>)});
}


template<class T>
T deserializeTrivial(const char*& pc) {
    T t;
    memcpy(&t, pc, sizeof(T));
    pc += sizeof(T);
    return t;
}

#endif //PHASEFIELD_SERIALIZE_H
