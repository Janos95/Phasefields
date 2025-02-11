//
// Created by janos on 11/2/20.
//

#pragma once

namespace Phasefield {

template<class T> struct RemoveReference      {typedef T type;};
template<class T> struct RemoveReference<T&>  {typedef T type;};
template<class T> struct RemoveReference<T&&> {typedef T type;};
template<class T> using RemoveReferenceT = typename RemoveReference<T>::type;

#define MOVE(x) static_cast<RemoveReferenceT<decltype(x)>&&>(x)
#define FWD(x) static_cast<decltype(x)&&>(x)

}