#include <iostream>
#include <CLI11.hpp>
#include <omp.h>
#include <vector>

//{}[]

int main(int argc, char** argv) {
    CLI::App app{"OpenMP example"};

    unsigned nworkers = 4;
    app.add_option("-n,--nworkers", nworkers, "Amount of workers");

    CLI11_PARSE(app, argc, argv);

    std::vector<uint64_t> idx;
    idx.reserve(nworkers);

    omp_set_dynamic(0);
    #pragma omp parallel shared(idx) num_threads(nworkers)
    {
        uint64_t me = omp_get_thread_num();
        #pragma omp critical
        {
            idx.push_back(me);
            std::cout << "Worker " << me << " done\n";
        }
    }

    std::sort(idx.begin(), idx.end());
    std::cout << "Sorted idx:\n";
    for (auto id : idx) {
        std::cout << "\t" << id << '\n';
    } 
    return 0;
}
