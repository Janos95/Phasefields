//
// Created by Janos Meny on 9/4/19.
//

#pragma once

#include <Corrade/Utility/StlForwardString.h>
#include <Corrade/Containers/Pointer.h>

class ScopedTimer
{
public:

    explicit ScopedTimer(std::string name, bool verbose = false);

    ~ScopedTimer();

    static void printStatistics();

private:

    struct Impl;
    Corrade::Containers::Pointer<Impl> m_impl;
};
