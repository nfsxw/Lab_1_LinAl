#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include "gauss_solver.h"
#include "lu_solver.h"
#include "matrix_utils.h"
#include "timer.h"

#ifdef _WIN32
#define CLEAR_CMD "cls"
#else
#define CLEAR_CMD "clear"
#endif

void clearConsole() {
    std::system(CLEAR_CMD);
}

void run_experiment_single_system() {
    std::cout << "\n=== Experiment 4.1: Single System Solution ===" << std::endl;
    std::vector<int> sizes = {100, 200, 500, 1000};
    
    for (int n : sizes) {
        std::cout << "\nMatrix Size N = " << n << std::endl;
        
        Matrix A = MatrixUtils::generate_random_matrix(n, -1.0, 1.0, 42);
        Vector b = MatrixUtils::generate_random_vector(n, -1.0, 1.0, 42);
        
        Timer timer;
        double t_gauss_no, t_gauss_part, t_lu_decomp, t_lu_solve;
        Vector x_res;

        try {
            timer.start();
            x_res = GaussSolver::solve_no_pivot(A, b);
            t_gauss_no = timer.stop();
            double res_no = MatrixUtils::norm(MatrixUtils::subtract(MatrixUtils::multiply(A, x_res), b));

            timer.start();
            x_res = GaussSolver::solve_partial_pivot(A, b);
            t_gauss_part = timer.stop();
            double res_part = MatrixUtils::norm(MatrixUtils::subtract(MatrixUtils::multiply(A, x_res), b));

            timer.start();
            LUResult lu = LuSolver::decompose(A);
            t_lu_decomp = timer.stop();
            
            timer.start();
            x_res = LuSolver::solve(lu, b);
            t_lu_solve = timer.stop();
            double res_lu = MatrixUtils::norm(MatrixUtils::subtract(MatrixUtils::multiply(A, x_res), b));

            MatrixUtils::print_table_header();
            MatrixUtils::print_experiment_result("Gauss No Pivot", t_gauss_no, res_no);
            MatrixUtils::print_experiment_result("Gauss Partial Pivot", t_gauss_part, res_part);
            MatrixUtils::print_experiment_result("LU Decomposition", t_lu_decomp, 0.0);
            MatrixUtils::print_experiment_result("LU Substitution", t_lu_solve, 0.0);
            MatrixUtils::print_experiment_result("LU Total", t_lu_decomp + t_lu_solve, res_lu);

        } catch (const std::exception& e) {
            std::cerr << "Error during calculation for N=" << n << ": " << e.what() << std::endl;
        }
    }
}

void run_experiment_multiple_rhs() {
    std::cout << "\n=== Experiment 4.2: Multiple Right-Hand Sides ===" << std::endl;
    int n = 500;
    std::vector<int> k_values = {1, 10, 100};
    
    Matrix A = MatrixUtils::generate_random_matrix(n, -1.0, 1.0, 42);
    
    std::cout << std::left << std::setw(10) << "K" 
              << std::right << std::setw(20) << "Gauss Time (ms)" 
              << std::setw(20) << "LU Time (ms)" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    for (int k : k_values) {
        Timer timer;
        
        timer.start();
        for (int i = 0; i < k; ++i) {
            Vector b = MatrixUtils::generate_random_vector(n, -1.0, 1.0, i);
            GaussSolver::solve_partial_pivot(A, b);
        }
        double t_gauss = timer.stop();

        timer.start();
        LUResult lu = LuSolver::decompose(A);
        double t_decomp = timer.stop();
        
        timer.start();
        for (int i = 0; i < k; ++i) {
            Vector b = MatrixUtils::generate_random_vector(n, -1.0, 1.0, i);
            LuSolver::solve(lu, b);
        }
        double t_solve = timer.stop();
        double t_lu_total = t_decomp + t_solve;

        std::cout << std::left << std::setw(10) << k 
                  << std::right << std::fixed << std::setprecision(4) 
                  << std::setw(20) << t_gauss 
                  << std::setw(20) << t_lu_total << std::endl;
    }
}

void run_experiment_hilbert() {
    std::cout << "\n=== Experiment 4.3: Hilbert Matrix Accuracy ===" << std::endl;
    std::vector<int> sizes = {5, 10, 15};

    std::cout << std::left << std::setw(10) << "N" 
              << std::setw(25) << "Method" 
              << std::right << std::setw(15) << "Rel Error" 
              << std::setw(15) << "Residual" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (int n : sizes) {
        Matrix H = MatrixUtils::generate_hilbert_matrix(n);
        Vector x_exact(n, 1.0);
        Vector b = MatrixUtils::multiply(H, x_exact);

        auto check_method = [&](const std::string& name, const Vector& x_approx) {
            double rel_error = MatrixUtils::norm(MatrixUtils::subtract(x_approx, x_exact)) / MatrixUtils::norm(x_exact);
            double residual = MatrixUtils::norm(MatrixUtils::subtract(MatrixUtils::multiply(H, x_approx), b));
            std::cout << std::left << std::setw(10) << n 
                      << std::setw(25) << name 
                      << std::right << std::scientific << std::setprecision(4) 
                      << std::setw(15) << rel_error 
                      << std::setw(15) << residual << std::endl;
        };

        try {
            Vector x_no_pivot = GaussSolver::solve_no_pivot(H, b);
            check_method("Gauss No Pivot", x_no_pivot);
        } catch (...) {}

        try {
            Vector x_part_pivot = GaussSolver::solve_partial_pivot(H, b);
            check_method("Gauss Partial Pivot", x_part_pivot);
        } catch (...) {}

        try {
            LUResult lu = LuSolver::decompose(H);
            Vector x_lu = LuSolver::solve(lu, b);
            check_method("LU Decomposition", x_lu);
        } catch (...) {}
        
        std::cout << std::string(65, '-') << std::endl;
    }
}

void print_menu() {
    clearConsole();
    std::cout << "\n========================================" << std::endl;
    std::cout << "   SLAE Solver Comparison Menu" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "1. Exp 4.1: Single System Time Comparison" << std::endl;
    std::cout << "2. Exp 4.2: Multiple RHS Efficiency" << std::endl;
    std::cout << "3. Exp 4.3: Hilbert Matrix Accuracy" << std::endl;
    std::cout << "4. Run All Experiments" << std::endl;
    std::cout << "0. Exit" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Select option: ";
}

int main() {
    clearConsole();
    int choice;
    while (true) {
        print_menu();
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        switch (choice) {
            case 1: run_experiment_single_system(); break;
            case 2: run_experiment_multiple_rhs(); break;
            case 3: run_experiment_hilbert(); break;
            case 4:
                run_experiment_single_system();
                run_experiment_multiple_rhs();
                run_experiment_hilbert();
                break;
            case 0:
                clearConsole();
                std::cout << "Exiting..." << std::endl;
                return 0;
            default:
                std::cout << "Invalid option. Please try again." << std::endl;
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        std::cout << "\nPress Enter to return to menu...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin.get();
        clearConsole();
    }
}