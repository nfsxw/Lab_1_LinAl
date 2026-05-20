#include "matrix_utils.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace MatrixUtils {

    Matrix generate_random_matrix(int n, double min_val, double max_val, unsigned int seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(min_val, max_val);
        Matrix A(n, Vector(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = dis(gen);
            }
        }
        return A;
    }

    Vector generate_random_vector(int n, double min_val, double max_val, unsigned int seed) {
        std::mt19937 gen(seed + 1);
        std::uniform_real_distribution<double> dis(min_val, max_val);
        Vector b(n);
        for (int i = 0; i < n; ++i) {
            b[i] = dis(gen);
        }
        return b;
    }

    Matrix generate_hilbert_matrix(int n) {
        Matrix H(n, Vector(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                H[i][j] = 1.0 / (static_cast<double>(i + j + 1));
            }
        }
        return H;
    }

    Vector multiply(const Matrix& A, const Vector& x) {
        int n = static_cast<int>(A.size());
        Vector res(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res[i] += A[i][j] * x[j];
            }
        }
        return res;
    }

    double norm(const Vector& v) {
        double sum = 0.0;
        for (double val : v) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }

    Vector subtract(const Vector& a, const Vector& b) {
        int n = static_cast<int>(a.size());
        Vector res(n);
        for (int i = 0; i < n; ++i) {
            res[i] = a[i] - b[i];
        }
        return res;
    }

    void print_table_header() {
        std::cout << std::left << std::setw(30) << "Method" 
                  << std::right << std::setw(15) << "Time (ms)" 
                  << std::setw(15) << "Residual" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
    }

    void print_experiment_result(const std::string& name, double time_ms, double residual) {
        std::cout << std::left << std::setw(30) << name 
                  << std::right << std::fixed << std::setprecision(6) << std::setw(15) << time_ms 
                  << std::setw(15) << residual << std::endl;
    }
}