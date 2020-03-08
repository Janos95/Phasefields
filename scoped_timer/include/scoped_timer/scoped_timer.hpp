//
// Created by Janos Meny on 9/4/19.
//

#pragma once

#include <string>
#include <memory>

class ScopedTimer
{
public:

    explicit ScopedTimer(std::string name, bool verbose = false);

    ~ScopedTimer();

    static void printStatistics();

private:

    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
