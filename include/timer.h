#pragma once
#include <chrono>

class Timer {
public:
    void start();
    double stop();
private:
    std::chrono::high_resolution_clock::time_point start_time;
};