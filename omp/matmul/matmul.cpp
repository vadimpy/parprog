#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <chrono>

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
        
        /*
        #pragma omp critical
        std::cout << "Thread " << me << ": range [" << begin << "; " << end
                  << ")\n";
        */

        for (size_t i = begin; i < end; ++i) {
            size_t row = i / k;
            size_t col = i % k;
            res[row][col] = inner(A[row], B[col]);
        }
    }
    return res;
}

std::vector<std::vector<double>> a_block_matmul(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b, int threads_num) {
    const size_t block_size = 16;
    size_t n = a.size();
    size_t m = a[0].size();
    size_t k = b.size();

    std::vector<std::vector<double>> res(n);
    std::for_each(res.begin(), res.end(), [k](std::vector<double>& row) {row.resize(k);});
    for (auto &row : res)
        for (auto &e : row)
            e = 0.0;

    #pragma omp parallel for schedule(static, 1)
    for (size_t block_i = 0; block_i < n; block_i += block_size)
        for (size_t block_j = 0; block_j < m; block_j += block_size)
            for (size_t block_k = 0; block_k < k; block_k += block_size) {
                size_t i_end = std::min(block_i + block_size, n);
                size_t j_end = std::min(block_j + block_size, m);
                size_t k_end = std::min(block_k + block_size, k);
                for (size_t i = block_i; i < i_end; ++i)
                    for (size_t j = block_j; j < j_end; ++j)
                        for (size_t k = block_k; k < k_end; ++k)
                            res[i][k] += a[i][j] * b[k][j];
            }

    return res;
}

std::vector<std::vector<double>> block_matmul(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b, int threads_num) {
    const size_t block_size = 16;
    size_t n = a.size();
    size_t m = a[0].size();
    size_t k = b.size();

    std::vector<std::vector<double>> res(n);
    std::for_each(res.begin(), res.end(), [k](std::vector<double>& row) {row.resize(k);});
    for (auto &row : res)
        for (auto &e : row)
            e = 0.0;

    #pragma omp parallel for schedule(static, 1)
    for (size_t i = 0; i < n; ++i)
        for (size_t block_j = 0; block_j < m; block_j += block_size)
            for (size_t block_k = 0; block_k < k; block_k += block_size) {
                size_t j_end = std::min(block_j + block_size, m);
                size_t k_end = std::min(block_k + block_size, k);
                for (size_t j = block_j; j < j_end; ++j)
                    for (size_t k = block_k; k < k_end; ++k)
                        res[i][k] += a[i][j] * b[k][j];
            }
    return res;
}

void print_matrix(std::vector<std::vector<double>> &m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m[0].size(); ++j)
            std::cout << m[i][j] << ' ';
        std::cout << '\n';
    }
}

void init_matrix(std::vector<std::vector<double>> &m) {
    for (auto &row : m) {
        std::transform(row.begin(), row.end(), row.begin(), [](double x) {
            return static_cast<double>(std::rand() % 1000) / 1000.0;
        });
    }
}

bool mat_cmp(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b) {
    if (a.size() != b.size() || a[0].size() != b[0].size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j)
            if (a[i][j] != b[i][j])
                return false;
    return true;
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
    // std::cout << "A:\n";
    // print_matrix(A);

    std::cout << '\n';

    std::vector<std::vector<double>> B(k);
    for (auto& row : B)
        row.resize(m);
    init_matrix(B);
    // std::cout << "B:\n";
    // print_matrix(B);

    std::chrono::high_resolution_clock c;
    auto t1 = c.now();
    auto res = block_mode ? block_matmul(A, B, num_threads) : matmul(A, B, num_threads);
    auto t2 = c.now();
    auto res1 = a_block_matmul(A, B, num_threads);
    auto t3 = c.now();

    std::cout << "5 loops mode time " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms\n";
    std::cout << "6 loops mode time " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << " ms\n";
    // std::cout << "Res:\n";
    // print_matrix(res);
    t1 = c.now();
    auto true_res = matmul(A, B, 1);
    std::cout << "validation time " << std::chrono::duration_cast<std::chrono::milliseconds>(c.now() - t1).count() << " ms\n";
    std::cout << "res: " << mat_cmp(res, true_res) << std::endl;

    return 0;
}
