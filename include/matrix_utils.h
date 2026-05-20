#pragma once
#include <vector>
#include <string>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

namespace MatrixUtils {
    Matrix generate_random_matrix(int n, double min_val, double max_val, unsigned int seed);
    Vector generate_random_vector(int n, double min_val, double max_val, unsigned int seed);
    Matrix generate_hilbert_matrix(int n);
    Vector multiply(const Matrix& A, const Vector& x);
    double norm(const Vector& v);
    Vector subtract(const Vector& a, const Vector& b);
    void print_table_header();
    void print_experiment_result(const std::string& name, double time_ms, double residual = 0.0);
}