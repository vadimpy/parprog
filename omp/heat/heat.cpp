#include <CLI11.hpp>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>

#define chunk(jobs, size) size / jobs + (size % jobs != 0)
#define begin(number, chunk) number * chunk
#define end(begin, chunk, size) std::min(begin + chunk, size)
#define correct_begin(index) if (index != 0) index += 1
#define correct_end(index, size) if (index != size) index -= 1

// u[i][j] = alpha * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j]) + u[i][j]

void heat_step_with_mem(std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &u_res, size_t jobs, double alpha) {

    #define d2u_dx2(u, i, j) (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])

    size_t n = u.size();
    size_t m = u[0].size();

    #pragma omp parallel for num_threads(jobs) collapse(2)
    for (size_t i = 1; i < n - 1; ++i)
        for (size_t j = 1; j < m - 1; ++j)
                u_res[i][j] = alpha * d2u_dx2(u, i, j) + u[i][j];

    #undef d2u_dx2
}

inline void process_first_row(
    std::vector<std::vector<double>> &u,
    std::vector<double> &prev_row,
    double alpha,
    size_t m
) {
    for (size_t j = 0; j < m; ++j)
        prev_row[j] = u[0][j];

    double first = u[0][0];
    u[0][0] = 2 * alpha * (u[1][0] + u[0][1] - 2 * first) + first;

    double prev = first;
    for (size_t j = 1; j < m - 1; ++j) {
        double cur = u[0][j];
        u[0][j] = 4 * alpha / 3 * (prev + u[1][j] + u[0][j+1] - 3 * cur) + cur;
        prev = cur;
    }

    double last = u[0][m-1];
    u[0][m-1] = 2 * alpha * (u[1][m-1] + prev - 2 * last) + last;
}

inline void init_prev_row(
    std::vector<std::vector<double>> &u,
    std::vector<double> &prev_row,
    std::vector<std::vector<double>> &edge_rows_buf,
    size_t m,
    size_t t_num
) {
    size_t prev_row_buf_idx = (t_num - 1) * 2;
    for (size_t j = 0; j < m; ++j)
        prev_row[j] = edge_rows_buf[prev_row_buf_idx][j];
}

inline void process_mid_rows(
    std::vector<std::vector<double>> &u,
    std::vector<double> &prev_row,
    size_t begin,
    size_t end,
    size_t m,
    double alpha
) {
    for (size_t i = begin; i < end - 1; ++i) {

        double first = u[i][0];
        u[i][0] = 4 * alpha / 3 * (prev_row[0] + u[i][1] + u[i+1][0] - 3 * first) + first;
        double prev = first;

        for (size_t j = 1; j < m - 1; ++j) {
            double cur = u[i][j];
            u[i][j] = alpha * (prev + prev_row[j] + u[i+1][j] + u[i][j+1] - 4 * cur) + cur;
            prev = prev_row[j] = cur;
        }

        double last = u[i][m-1];
        u[i][m-1] = 4 * alpha / 3 * (prev + prev_row[m-1] + u[i+1][m-1] - 3 * last) + last;
    }
}

inline void process_very_last_row(
    std::vector<std::vector<double>> &u,
    std::vector<double> &prev_row,
    size_t n,
    size_t m,
    double alpha
) {
    double first = u[n-1][0];
    u[n-1][0] = 2 * alpha * (prev_row[0] + u[n-1][1] - 2 * first) + first;

    double prev = first;
    for (size_t j = 1; j < m - 1; ++j) {
        double cur = u[n-1][j];
        u[n-1][j] = 4 * alpha / 3 * (prev + prev_row[j] + u[n-1][j+1] - 3 * cur) + cur;
        prev = cur;
    }
    double last = u[n-1][m-1];
    u[n-1][m-1] = 2 * alpha * (prev + prev_row[m-1] - 2 * last) + last;
}

inline void process_chunk_last_row(
    std::vector<std::vector<double>> &u,
    std::vector<double> &prev_row,
    std::vector<std::vector<double>> edge_rows_buf,
    double alpha,
    size_t end,
    size_t m,
    size_t t_num
) {
    std::vector<double>& next_row = edge_rows_buf[t_num * 2 + 1];

    double first = u[end-1][0];
    u[end-1][0] = 4 * alpha / 3 * (prev_row[0] + next_row[0] + u[end-1][1] - 3 * first) + first;

    double prev = first;
    for (size_t j = 1; j < m - 1; ++j) {
        double cur = u[end-1][j];
        u[end-1][j] = alpha * (prev + prev_row[j] + next_row[j] + u[end-1][j+1] - 4 * cur) + cur;
        prev = cur;
    }

    double last = u[end-1][m-1];
    u[end-1][m-1] = 4 * alpha / 3 * (prev + prev_row[m-1] + next_row[m-1] - 3 * last) + last;

}

void heat_step(std::vector<std::vector<double>> &u, size_t jobs, double alpha) {
    size_t n = u.size();
    size_t m = u[0].size();

    size_t chunk = chunk(jobs, n);
    size_t pairs_amount = jobs - 1;

    std::vector<std::vector<double>> edge_rows_buf(2 * pairs_amount);
    for (auto &x: edge_rows_buf)
        x.resize(m);

    for (size_t i = 0, k = chunk - 1; i < pairs_amount; i += 1, k += chunk)
        for (size_t j = 0; j < m; ++j) { 
            edge_rows_buf[2*i][j] = u[k][j];
            edge_rows_buf[2*i+1][j] = u[k+1][j];
        }

    #pragma omp parallel num_threads(jobs) shared(u, edge_rows_buf)
    {
        auto t_num = omp_get_thread_num();
        std::vector<double> prev_row(m);

        size_t begin = begin(t_num, chunk);
        size_t end = end(begin, chunk, n);

        if (begin == 0) {
            process_first_row(u, prev_row, alpha, m);
            ++begin;
        } else
            init_prev_row(u, prev_row, edge_rows_buf, m, t_num);

        process_mid_rows(u, prev_row, begin, end, m, alpha);

        if (end == n)
            process_very_last_row(u, prev_row, n, m, alpha);
        else
            process_chunk_last_row(u, prev_row, edge_rows_buf, alpha, end, m, t_num);
    }
    
}

int main(int argc, char **argv) {

    size_t jobs = 2;
    size_t m = 100;
    size_t n = 100;
    size_t frames = 100;
    double alpha = 0.24;

    CLI::App app("Parallel numerical solution for heat equation");
    app.add_option("-j, --jobs", jobs, "Number of threads (defaults to 2)");
    app.add_option("-m", m, "M size (defaults to 100)");
    app.add_option("-n", n, "N size (defaults to 100)");
    app.add_option("-a, --alpha", alpha, "Heat coefficient (defaults to 0.24, please, don't set if greater than 0.25 :) )");
    app.add_option("--frames", frames, "Amount of frames (defaults to 100)");
    CLI11_PARSE(app, argc, argv);

    std::vector<std::vector<double>> u(n);
    for(size_t i = 0; i < n; ++i)
        u[i].resize(m);

    for (size_t j = 0; j < m; ++j)
        u[n-1][j] = 400;
    for (size_t i = 0; i < n; ++i) {
        u[i][0] = 400;
        u[i][m-1] = 400;
    }

    for (size_t i = n/2 - n/8; i < n/2 + n/8; ++i)
        u[i][m/2] = 400;

    for (size_t j = m/2 - m/8; j < m/2 + m/8; ++j)
        u[n/2][j] = 400;

    for (size_t i = 0; i < frames; ++i) {
        heat_step(u, jobs, alpha);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j)
                std::cout << u[i][j] << ' ';
            std::cout << '\n';
        }
        std::cout << '\n';
    }
    return 0;
}
