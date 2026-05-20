#pragma once
#include "matrix_utils.h"

namespace GaussSolver {
    Vector solve_no_pivot(const Matrix& A, const Vector& b);
    Vector solve_partial_pivot(const Matrix& A, const Vector& b);
}