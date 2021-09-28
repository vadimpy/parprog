#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>

typedef std::vector<std::vector<double>> matrix;

double inner(std::vector<double>& a, std::vector<double>& b) {
    double res = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        res += a[i] * b[i];
    return res;
}

std::vector<std::vector<double>> matmul(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B, int threads_num) {
    size_t n = A.size();
    size_t m = A[0].size();
    size_t k = B.size();
    size_t len = n*k;
    size_t chunk_size = len / threads_num + (len % threads_num != 0);

    std::vector<std::vector<double>> res(n);
    std::for_each(res.begin(), res.end(), [k](std::vector<double>& row) {row.resize(k);});

    #pragma omp parallel num_threads(threads_num)
    {
        int me = omp_get_thread_num();
        size_t begin = me * chunk_size;
        size_t end = std::min(begin + chunk_size, len);
        
        #pragma omp critical
        std::cout << "Thread " << me << ": range [" << begin << "; " << end
                  << ")\n";

        for (size_t i = begin; i < end; ++i) {
            size_t row = i / k;
            size_t col = i % k;
            res[row][col] = inner(A[row], B[col]);
        }
    }
    return res;
}

std::vector<std::vector<double>> block_matmul(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B, int threads_amount) {
    size_t len = A.size() * A.size();

}

void print_matrix(std::vector<std::vector<double>> &m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m[0].size(); ++j)
            std::cout << m[i][j] << ' ';
        std::cout << '\n';
    }
}

void init_matrix(std::vector<std::vector<double> > &m) {
    for (auto &row : m) {
        std::transform(row.begin(), row.end(), row.begin(), [](double x) {
            return static_cast<double>(std::rand() % 1000) / 1000.0;
        });
    }
}

int main(int argc, char **argv) {

    size_t n = 256;
    size_t m = 256;
    size_t k = 256;
    bool block_mode = false;
    size_t num_threads = 1;
    unsigned seed = 0;

    CLI::App app{"Matrix multiplication with 2 modes: simple and block (NxM * MxK = NxK)"};
    app.add_option("-n", n,
                   "N size");
    app.add_option("-m", m,
                   "M size");
    app.add_option("-k", k,
                   "K size");
    app.add_flag(
        "--block", block_mode,
        "Enables block multiplication");
    app.add_option("-j, --jobs", num_threads,
                   "Number of threads (defaults to 1)");
    app.add_option("--seed", seed, "Seed for random generator (defaults to 0)");
    CLI11_PARSE(app, argc, argv);

    std::srand(seed);

    std::vector<std::vector<double>> A(n);
    for (auto& row : A)
        row.resize(m);
    init_matrix(A);
    std::cout << "A:\n";
    print_matrix(A);

    std::cout << '\n';

    std::vector<std::vector<double>> B(k);
    for (auto& row : B)
        row.resize(m);
    init_matrix(B);
    std::cout << "B:\n";
    print_matrix(B);

    auto res = block_mode ? block_matmul(A, B, num_threads) : matmul(A, B, num_threads);
    std::cout << "Res:\n";
    print_matrix(res);

    return 0;
}
