#pragma once
#include <Eigen/Dense>
#include <ostream>
#include <vector>

enum class SolveStatus { OPTIMAL = 1, UNBOUNDED = -1, INFEASIBLE = -2 };

inline std::ostream& operator<<(std::ostream& os, SolveStatus s) {
    switch (s) {
        case SolveStatus::OPTIMAL:    return os << "OPTIMAL";
        case SolveStatus::UNBOUNDED:  return os << "UNBOUNDED";
        case SolveStatus::INFEASIBLE: return os << "INFEASIBLE";
    }
    return os << static_cast<int>(s);
}

struct LinearProgram {
    Eigen::MatrixXd A;
    Eigen::VectorXd b, c;

    LinearProgram(const Eigen::MatrixXd& A_,
                  const Eigen::VectorXd& b_,
                  const Eigen::VectorXd& c_)
        : A(A_), b(b_), c(c_) {}
};

struct BasisState {
    std::vector<int> basic, nonbasic;
};

struct SimplexResult {
    Eigen::VectorXd x;
    double          z;
    SolveStatus     status;
    int             iterations;
};
