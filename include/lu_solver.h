#pragma once
#include "matrix_utils.h"

struct LUResult {
    Matrix L;
    Matrix U;
};

namespace LuSolver {
    LUResult decompose(const Matrix& A);
    Vector forward_substitution(const Matrix& L, const Vector& b);
    Vector backward_substitution(const Matrix& U, const Vector& y);
    Vector solve(const LUResult& lu, const Vector& b);
}