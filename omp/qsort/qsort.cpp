#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <chrono>

enum InputMode {
    RAND = 0,
    STDIN = 1,
    FILE_x = 2
};

void qsort(std::vector<int64_t> &a, size_t begin, size_t end, uint16_t current_depth, uint16_t max_depth) {
    if (begin >= end) return;
    size_t len = end + 1 - begin;
    if (len == 2) {
        if (a[begin] > a[end])
            std::swap(a[begin], a[end]);
        return;
    }
    int64_t pivot = a[std::rand() % len + begin];
    int64_t pivot_index = begin;
    for (size_t i = begin; i <= end; ++i)
        if (a[i] < pivot) {
            std::swap(a[pivot_index], a[i]);
            ++pivot_index;
        }
    if (current_depth == max_depth) {
        if (pivot_index > begin)
            qsort(a, begin, pivot_index - 1, current_depth, max_depth);
        if (a[pivot_index] == pivot)
            ++pivot_index;
        qsort(a, pivot_index, end, current_depth, max_depth);
    } else
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                if (pivot_index > begin)
                    qsort(a, begin, pivot_index - 1, current_depth + 1, max_depth);
            }
            #pragma omp section
            {
                if (a[pivot_index] == pivot)
                    ++pivot_index;
                qsort(a, pivot_index, end, current_depth + 1, max_depth);
            }
        }
}

bool check_sorted(const std::vector<int64_t> &a) {
    for(size_t i = 1; i < a.size(); ++i)
        if (a[i-1] > a[i])
            return false;
    return true;
}

int main(int argc, char **argv) {
    size_t n = 256;
    size_t jobs = 1;
    uint64_t seed = 0;
    int64_t a = -100;
    int64_t b = 100;
    uint16_t depth = 4;
    bool verbose = false;
    InputMode in_mode = InputMode::RAND;

    CLI::App app{"Parallel qsort"};
    app.add_option("-n", n, "length of array to be generated (defaults to 256)");
    app.add_option("--seed", seed, "Seed (defaults to 0)");
    app.add_option("-a", a, "Left edge of random range (defaults to -100)");
    app.add_option("-b", b, "Right edge of random range (defaults to 100)");
    app.add_option("--depth", depth, "Max depth for thread creation");
    app.add_option("--input-mode", in_mode, "Input mode (0 - random, 1 - stdin, 2 - input.txt)");
    app.add_flag("--verbose", verbose, "Output sorted array");

    CLI11_PARSE(app, argc, argv);

    std::srand(seed);

    omp_set_max_active_levels(depth);

    std::vector<int64_t> arr;
    switch (in_mode) {
        case InputMode::RAND:
            arr.resize(n);
            for (size_t i = 0; i < n; ++i)
                arr[i] = (std::rand() % (b - a + 1) + a);
        break;
        case InputMode::STDIN:
            arr.resize(n);
            for (size_t i = 0; i < n; ++i)
                std::cin >> arr[i];
        break;
        case InputMode::FILE_x:
            std::ifstream f("input.txt");
            while (!f.eof()) {
                int64_t x;
                f >> x;
                arr.push_back(x);
            }
        break;
    }

    qsort(arr, 0, arr.size() - 1, 0, depth);
    std::cout << (check_sorted(arr) ? "array is sorted\n" : "array is not sorted\n");
    if (verbose) {
        for (auto a: arr)
            std::cout << a << ' ';
        std::cout << '\n';
    }
    return 0;
}
