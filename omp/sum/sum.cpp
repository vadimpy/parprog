#include <iostream>
#include <CLI11.hpp>
#include <omp.h>
#include <vector>

//{}[]

int main(int argc, char** argv) {
    unsigned n = 10;

    CLI::App app{"OpenMP sum from 1 to M"};
    app.add_option("-n", n, "Upper edge of sum range");
    CLI11_PARSE(app, argc, argv);

    unsigned sum = 0;

    #pragma omp parallel for reduction(+:sum) schedule(auto)
    for (unsigned i = 1; i <= n; ++i)
        sum += i;

    std::cout << sum << std::endl;
}
