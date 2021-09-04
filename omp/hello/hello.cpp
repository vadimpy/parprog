#include <iostream>
#include <CLI11.hpp>
#include <omp.h>
#include <vector>
#include <string>

//{}[]

const char hello_world[] = "hello world";

int main(int argc, char** argv) {
    CLI::App app{"OpenMP hello world"};

    const auto nworkers = sizeof(hello_world);

    unsigned next_id = 0;
    CLI11_PARSE(app, argc, argv);

    #pragma omp parallel shared(hello_world, next_id) num_threads(nworkers)
    {
        auto me = omp_get_thread_num();
        while (next_id != me);
        std::cout << hello_world[me];
        ++next_id;
    }
    std::cout << std::endl;
    return 0;
}
