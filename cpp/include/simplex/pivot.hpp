#pragma once
#include <optional>
#include <Eigen/Dense>

struct PivotRule {
    virtual std::optional<int> entering_index(const Eigen::VectorXd& cr,
                                               double tol = 1e-10) const = 0;
    virtual ~PivotRule() = default;
};

struct LargestCoefficient : PivotRule {
    std::optional<int> entering_index(const Eigen::VectorXd& cr,
                                      double tol = 1e-10) const override;
};

struct Bland : PivotRule {
    std::optional<int> entering_index(const Eigen::VectorXd& cr,
                                      double tol = 1e-10) const override;
};

std::optional<int> leaving_index(const Eigen::VectorXd& xb,
                                  const Eigen::VectorXd& d,
                                  const Eigen::MatrixXd& BinvN,
                                  double tol = 1e-10);
