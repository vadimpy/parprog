#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <chrono>

bool is_simple(uint64_t num) {
    if (num <= 2)
        return true;

    for (size_t i = 2; i * i < num; ++i)
        if (num % i == 0)
            return false;
    return true;
}

std::vector<uint64_t> simple_search(size_t n, size_t jobs) {
    size_t chunk = n / jobs + (n % jobs != 0);
    std::vector<uint64_t> res;
    #pragma omp parallel for num_threads(jobs)
    for (size_t i = 1; i <= n; ++i)
        if (is_simple(i))
        #pragma omp critical
            res.push_back(i);
    return res;
}

std::vector<bool> sieve_search(size_t n, size_t jobs) {
    std::vector<bool> v(n);
    for(size_t i = 0; i < n; ++i)
        v[i] = true;

    v[0] = false;
    v[1] = false;
    v[2] = true;

    size_t start = 2;

    while (start < n) {
        for(;!v[start];++start);
        #pragma omp parallel for num_threads(jobs)
        for (size_t i = 2 * start; i < n; i += start)
            v[i] = false;
        ++start;
    }

    return v;
}

template<typename T>
void print_vec(const std::vector<T>& v) {
    for (auto &x: v)
        std::cout << x << ' ';
    std::cout << std::endl;
}

void filtered_print(const std::vector<bool>& v) {
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i])
            std::cout << i << ' ';
    std::cout << std::endl;
}

int main(int argc, char **argv) {

    size_t n = 256;
    bool sieve_mode = false;
    size_t jobs = 1;
    bool silent = false;

    CLI::App app{"Simple numbers from 1 to N"};
    app.add_option("-n", n, "N");
    app.add_flag(
        "--sieve", sieve_mode,
        "Enables block multiplication");
    app.add_option("-j, --jobs", jobs,
                   "Number of threads (defaults to 1)");
    app.add_flag("--silent", silent, "No stdout");
    CLI11_PARSE(app, argc, argv);

    if (sieve_mode) {
        auto res = sieve_search(n, jobs);
        if (!silent)
            filtered_print(res);
    } else {
        auto res = simple_search(n, jobs);
        if (!silent)
            print_vec(res);
    }
    
    return 0;
}
