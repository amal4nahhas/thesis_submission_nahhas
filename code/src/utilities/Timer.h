#pragma once

#include <ctime>
#include <chrono>

class Timer
{
   public:
    using ClockType = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<ClockType>;

    Timer() noexcept : start_(ClockType::now()) {}

    void reset() noexcept { start_ = ClockType::now(); }

    double elapsedSeconds() const noexcept 
    { 
        return std::chrono::duration<double>(ClockType::now() - start_).count(); 
    }

   private:
    TimePoint start_;
};



class TimerCtime
{
   public:
    using ClockType = std::clock_t;
    static constexpr double TicksPerSecond = static_cast<double>(CLOCKS_PER_SEC);

    TimerCtime() noexcept : start_(std::clock()) {}

    void reset() noexcept { start_ = std::clock(); }

    ClockType elapsedTicks() const noexcept { return std::clock() - start_; }

    double elapsedSeconds() const noexcept { return elapsedTicks() / TicksPerSecond; }

    static double toSeconds(ClockType ticks) noexcept { return ticks / TicksPerSecond; }

   private:
    ClockType start_;
};

