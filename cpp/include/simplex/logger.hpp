#pragma once
#include <ostream>
#include <string>
#include "types.hpp"

enum class LogLevel { SILENT, INFO, DEBUG };

class SimplexLogger {
public:
    explicit SimplexLogger(LogLevel level = LogLevel::SILENT,
                           std::ostream& stream = nullstream());

    void log_problem(const LinearProgram& lp) const;
    void log_phase(int phase) const;
    void log_iteration(int it, const BasisState& basis,
                       const Eigen::VectorXd& xb, double z) const;
    void log_result(const SimplexResult& result) const;

private:
    LogLevel       level_;
    std::ostream&  out_;

    void write(const std::string& s) const;

    static std::ostream& nullstream();
};
