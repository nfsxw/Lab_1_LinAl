#include "lu_solver.h"
#include <stdexcept>

namespace LuSolver {

    LUResult decompose(const Matrix& A) {
        int n = static_cast<int>(A.size());
        Matrix L(n, Vector(n, 0.0));
        Matrix U = A; 

        for (int i = 0; i < n; ++i) {
            L[i][i] = 1.0;
        }

        for (int k = 0; k < n - 1; ++k) {
            if (std::abs(U[k][k]) < 1e-12) throw std::runtime_error("Zero pivot in LU");
            for (int i = k + 1; i < n; ++i) {
                double factor = U[i][k] / U[k][k];
                L[i][k] = factor;
                for (int j = k; j < n; ++j) {
                    U[i][j] -= factor * U[k][j];
                }
            }
        }
        return {L, U};
    }

    Vector forward_substitution(const Matrix& L, const Vector& b) {
        int n = static_cast<int>(L.size());
        Vector y(n);
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int j = 0; j < i; ++j) {
                sum -= L[i][j] * y[j];
            }
            y[i] = sum / L[i][i];
        }
        return y;
    }

    Vector backward_substitution(const Matrix& U, const Vector& y) {
        int n = static_cast<int>(U.size());
        Vector x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[i];
            for (int j = i + 1; j < n; ++j) {
                sum -= U[i][j] * x[j];
            }
            x[i] = sum / U[i][i];
        }
        return x;
    }

    Vector solve(const LUResult& lu, const Vector& b) {
        Vector y = forward_substitution(lu.L, b);
        return backward_substitution(lu.U, y);
    }
}