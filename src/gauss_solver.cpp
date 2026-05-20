#include "gauss_solver.h"
#include <cmath>
#include <stdexcept>

namespace GaussSolver {

    namespace {
        Vector forward_elimination_no_pivot(Matrix A, Vector b) {
            int n = static_cast<int>(A.size());
            for (int k = 0; k < n - 1; ++k) {
                if (std::abs(A[k][k]) < 1e-12) throw std::runtime_error("Zero pivot");
                for (int i = k + 1; i < n; ++i) {
                    double factor = A[i][k] / A[k][k];
                    for (int j = k; j < n; ++j) {
                        A[i][j] -= factor * A[k][j];
                    }
                    b[i] -= factor * b[k];
                }
            }
            
            Vector x(n);
            for (int i = n - 1; i >= 0; --i) {
                double sum = b[i];
                for (int j = i + 1; j < n; ++j) {
                    sum -= A[i][j] * x[j];
                }
                x[i] = sum / A[i][i];
            }
            return x;
        }

        Vector forward_elimination_partial_pivot(Matrix A, Vector b) {
            int n = static_cast<int>(A.size());
            for (int k = 0; k < n - 1; ++k) {
                int max_row = k;
                double max_val = std::abs(A[k][k]);
                for (int i = k + 1; i < n; ++i) {
                    if (std::abs(A[i][k]) > max_val) {
                        max_val = std::abs(A[i][k]);
                        max_row = i;
                    }
                }
                
                if (max_row != k) {
                    std::swap(A[k], A[max_row]);
                    std::swap(b[k], b[max_row]);
                }

                if (std::abs(A[k][k]) < 1e-12) throw std::runtime_error("Zero pivot");

                for (int i = k + 1; i < n; ++i) {
                    double factor = A[i][k] / A[k][k];
                    for (int j = k; j < n; ++j) {
                        A[i][j] -= factor * A[k][j];
                    }
                    b[i] -= factor * b[k];
                }
            }

            Vector x(n);
            for (int i = n - 1; i >= 0; --i) {
                double sum = b[i];
                for (int j = i + 1; j < n; ++j) {
                    sum -= A[i][j] * x[j];
                }
                x[i] = sum / A[i][i];
            }
            return x;
        }
    }

    Vector solve_no_pivot(const Matrix& A, const Vector& b) {
        return forward_elimination_no_pivot(A, b);
    }

    Vector solve_partial_pivot(const Matrix& A, const Vector& b) {
        return forward_elimination_partial_pivot(A, b);
    }
}