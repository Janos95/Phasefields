//
// Created by janos on 19.10.19.
//

#include "scoped_timer.hpp"

#include <fmt/format.h>

#include <mutex>
#include <unordered_map>
#include <chrono>
#include <cmath>

using my_clock = std::chrono::steady_clock;
using my_duration = std::chrono::duration<long double, std::nano>;
using time_point_t = std::chrono::time_point<my_clock>;

template<class Ratio>
std::string toSI()
{
    if constexpr(std::is_same_v<std::nano, Ratio>)
        return "nano seconds";
    if constexpr(std::is_same_v<std::micro, Ratio>)
        return "micro seconds";
    if constexpr(std::is_same_v<std::milli, Ratio>)
        return "milli seconds";
    if constexpr(std::is_same_v<std::ratio<1>, Ratio>)
        return "seconds";
    if constexpr(std::is_same_v<std::ratio<60>, Ratio>)
        return "minutes";
    if constexpr(std::is_same_v<std::ratio<3600>, Ratio>)
        return "hours";
}

using Ratio = std::ratio<1>;
using user_dur = std::chrono::duration<long double, Ratio>;

struct TimingInfo
{
    my_duration mean{0};
    long double M2 = 0;
    std::size_t count = 0;
};

struct ScopedTimer::Impl{
    std::string name_;
    time_point_t start_;
    bool verbose_;
};

std::unordered_map<std::string, TimingInfo> log_;
std::mutex mutex_;

ScopedTimer::ScopedTimer(std::string name, bool verbose):
    m_impl(std::make_unique<Impl>(Impl{std::move(name), my_clock::now(), verbose}))
{
}

ScopedTimer::~ScopedTimer()
{
    const auto end = my_clock::now();
    my_duration time = end - m_impl->start_;
    if(m_impl->verbose_)
        fmt::print("{} took {} {}\n", m_impl->name_, user_dur{time}.count(), toSI<Ratio>());

    std::lock_guard l(mutex_);

    auto& [mean, M2, count] = log_[m_impl->name_];
    ++count;
    auto delta1 = time - mean;
    mean += delta1 / count;
    auto delta2 = time - mean;
    M2 += delta1.count() * delta2.count();
}

void ScopedTimer::printStatistics()
{

    std::lock_guard l(mutex_);
    for(const auto& [name, timingInfo]: log_)
    {
        const auto& [mean, M2, count] = timingInfo;
        user_dur standardDeviation = my_duration{std::sqrt(M2/(count - 1))};
        user_dur meanUser = mean;
        auto unit = toSI<Ratio>();
        fmt::print("{3}: Mean {0} {2}, Standard Deviation {1} {2}\n", meanUser.count(), standardDeviation.count(), unit, name);
    }
}
