//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_TIMER_H
#define MLCDEC_TIMER_H


#include "../Header.h"

class Timer {
public:
    Timer() {
        start = omp_get_wtime();
    }

    ~Timer() = default;

    void Stop() {
        end = omp_get_wtime();
    }

    [[nodiscard]] double GetTimeInSec() const {
        return end - start;
    }

private:
    double start{};
    double end{};
};

#endif //MLCDEC_TIMER_H
