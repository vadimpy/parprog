#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>

double comp_exp(double x0, size_t len, size_t jobs) {
    std::cout << "Exp" << '\n'; 
    size_t chunk = len / jobs + (len % jobs != 0);
    double res = 0.0;
    #pragma omp parallel num_threads(jobs)
    {
        auto me = omp_get_thread_num();
        size_t begin = me * chunk;
        size_t end = std::min(begin + chunk, len);

        double fact = 1;
        double x = 1;
        for (size_t i = 1; i <= begin; ++i) {
            fact /= i;
            x *= x0;
        }

        double tmp = x * fact;
        for (size_t i = begin + 1; i < end; ++i) {
            fact /= i;
            x *= x0;
            tmp += x * fact;
        }

        #pragma omp atomic
        res += tmp;
    }
    return res;
}

double comp_cos(double x0, size_t len, size_t jobs) {
    std::cout << "Cos" << '\n';
    size_t chunk = len / jobs + (len % jobs != 0);
    double res = 0.0;
    x0 = x0 * x0;
    #pragma omp parallel num_threads(jobs)
    {
        auto me = omp_get_thread_num();
        size_t begin = me * chunk;
        size_t end = std::min(begin + chunk, len);

        double fact = 1;
        double x = 1;
        for (size_t i = 1, j = 1; i <= begin; ++i, j += 2) {
            fact /= j*(j+1);
            x *= x0;
        }

        double sign = begin % 2 == 0 ? 1 : -1;
        double tmp = sign * x * fact;
        for (size_t i = begin + 1, j = i; i < end; ++i, j += 2) {
            fact /= j*(j+1);
            x *= x0;
            sign = -sign;
            tmp += sign * x * fact;
        }

        #pragma omp atomic
        res += tmp;
    }
    return res;
}

int main(int argc, char **argv) {

    uint16_t target = 0;
    double x = 0;
    size_t jobs = 4;
    size_t len = 16;

    CLI::App app("Series computation");
    app.add_option("-t, --target", target,
                   "Target to be computed ('--target 0' - exp(x), '--target 1' - cos(x))");
    app.add_option("-x", x, "x");
    app.add_option("-j, --jobs", jobs,
                   "Number of threads (defaults to 1)");
    app.add_option("-l, --len", len, "Length of series");
    CLI11_PARSE(app, argc, argv);

    if (target == 0) {
        double res = comp_exp(x, len, jobs);
        std::cout << "exp(" << x <<  ") ~ " << res << '\n';
    }
    else {
        double res = comp_cos(x, len, jobs);
        std::cout << "cos(" << x << ") ~ " << res << '\n';
    }

    return 0;
}
