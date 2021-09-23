#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>

void process_vec(std::vector<double> &data, const size_t threads_num) {

    size_t n = data.size();
    std::vector<double> tmp(n);

    #pragma omp parallel for num_threads(threads_num)
    for (size_t i = 1; i < n - 1; ++i) {
        tmp[i] = (data[i - 1] + data[i] + data[i + 1]) / 3.0;
    }

    tmp[0] = (data[0] + data[1]) / 2.0;
    tmp[n - 1] = (data[n - 1] + data[n - 2]) / 2.0;

    #pragma omp parallel for num_threads(threads_num)
    for (size_t i = 0; i < n; ++i) {
        data[i] = tmp[i];
    }
}

void process_vec_inplace(std::vector<double> &data, const size_t threads_num) {

    size_t n = data.size();
    size_t chunk_size = n / threads_num + (n % threads_num != 0);
    std::cout << "Chunk size: " << chunk_size << '\n';

    // save thread-dependent data in read-only place to prevent races in corners of chunks
    const std::vector<std::array<double, 2>> chunks_edges(threads_num + 1);
    for (int i = 1; i < threads_num; ++i) {
        size_t tmp = i * chunk_size;
        chunks_edges[i][0] = data[tmp - 1];
        chunks_edges[i][1] = data[tmp];
    }

    // trick: to deal with data[0] and data[n-1] regularly
    // use this equation: (x + y) / 2 == (x + y + a) / 3 => a = (x + y) / 2
    chunks_edges[0][0] = (data[0] + data[1]) / 2.0;
    chunks_edges[threads_num][1] = (data[n - 1] + data[n - 2]) / 2.0;

    #pragma omp parallel num_threads(threads_num)
    {
        int me = omp_get_thread_num();
        size_t begin = me * chunk_size;
        size_t end = std::min(begin + chunk_size, n);

        #pragma omp critical
        std::cout << "Thread " << me << ": chunk [" << begin << "; " << end
                  << ")\n";

        // save data[i] before changing to reuse it on the next iteration
        double prev = data[begin];
        data[begin] += chunks_edges[me][0] + data[begin+1];
        for (size_t i = begin + 1; i < end - 1; ++i) {
            double tmp = data[i];
            data[i] += prev + data[i+1];
            prev = tmp;
        }
        data[end-1] += chunks_edges[me+1][1] + prev;

        for (size_t i = begin; i < end; ++i)
            data[i] /= 3.0;
    }
}

int main(int argc, char **argv) {

    size_t n = 256;
    bool inplace = false;
    size_t num_threads = 1;
    unsigned seed = 0;

    CLI::App app{"Transformation a[i] = (a[i-1]+a[i]+a[i+1]) / 3 (inplace and extra-memory implementation)"};
    app.add_option("-s, --size", n,
                   "Length of vector to generate (defaults to 256)");
    app.add_flag(
        "--inplace", inplace,
        "Process vector inplace (e. g. without intermediate mem alloc)");
    app.add_option("-n, --nthreads", num_threads,
                   "Number of threads (defaults to 1)");
    app.add_option("--seed", seed, "Seed random generator (defaults to 0)");
    CLI11_PARSE(app, argc, argv);

    std::srand(seed);

    std::vector<double> data(n);
    std::transform(data.begin(), data.end(), data.begin(), [](double x) {
        return static_cast<double>(std::rand() % 1000) / 1000.0;
    });

    std::cout << "Seed: " << seed << "\nData size: " << data.size() << "\nGenerated data: ";
    for (auto x : data)
        std::cout << x << ' ';
    std::cout << "\n\n";

    if (inplace)
        process_vec_inplace(data, num_threads);
    else
        process_vec(data, num_threads);

    std::cout << "\nCalculated data: ";
    for (auto x : data)
        std::cout << x << ' ';
    std::cout << std::endl;

    return 0;
}
