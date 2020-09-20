//
// Created by Janos Meny on 9/4/19.
//

#pragma once

class ScopedTimer
{
public:

    explicit ScopedTimer(char const* name, bool verbose = false);

    ~ScopedTimer();

    static void printStatistics();

private:
    struct Impl;
    Impl* m_impl;
};
