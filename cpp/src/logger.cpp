#include "simplex/logger.hpp"
#include <iostream>
#include <sstream>
#include <string>

// A no-op streambuf used for the silent null stream
struct NullBuf : std::streambuf {
    int overflow(int c) override { return c; }
};

static NullBuf g_null_buf;
static std::ostream g_null_stream(&g_null_buf);

std::ostream& SimplexLogger::nullstream() { return g_null_stream; }

SimplexLogger::SimplexLogger(LogLevel level, std::ostream& stream)
    : level_(level), out_(stream) {}

void SimplexLogger::write(const std::string& s) const { out_ << s << '\n'; }

void SimplexLogger::log_problem(const LinearProgram& lp) const
{
    if (level_ == LogLevel::SILENT) return;
    write("=======================");
    write(" Revised Simplex Solver");
    write("=======================");
    {
        std::ostringstream ss;
        ss << "A =\n" << lp.A;
        write(ss.str());
    }
    {
        std::ostringstream ss;
        ss << "b = " << lp.b.transpose();
        write(ss.str());
    }
    {
        std::ostringstream ss;
        ss << "c = " << lp.c.transpose();
        write(ss.str());
    }
    write("");
}

void SimplexLogger::log_phase(int phase) const
{
    if (level_ == LogLevel::SILENT) return;
    std::string header = (phase == 1)
        ? "Phase 1 \xe2\x80\x94 Finding initial BFS"
        : "Phase 2 \xe2\x80\x94 Optimizing";
    write(header);
    write(std::string(header.size(), '-'));
}

void SimplexLogger::log_iteration(int it, const BasisState& basis,
                                   const Eigen::VectorXd& xb, double z) const
{
    if (level_ != LogLevel::DEBUG) return;
    std::ostringstream ss;
    ss << "  iter " << it << ":\n";
    ss << "    xb       = " << xb.transpose() << "\n";
    ss << "    basic    = [";
    for (std::size_t i = 0; i < basis.basic.size(); ++i)
        ss << basis.basic[i] << (i + 1 < basis.basic.size() ? ", " : "");
    ss << "]\n";
    ss << "    nonbasic = [";
    for (std::size_t i = 0; i < basis.nonbasic.size(); ++i)
        ss << basis.nonbasic[i] << (i + 1 < basis.nonbasic.size() ? ", " : "");
    ss << "]\n";
    ss << "    z        = " << z;
    write(ss.str());
}

void SimplexLogger::log_result(const SimplexResult& result) const
{
    if (level_ == LogLevel::SILENT) return;
    write("");
    std::string status_str;
    switch (result.status) {
        case SolveStatus::OPTIMAL:    status_str = "OPTIMAL";    break;
        case SolveStatus::UNBOUNDED:  status_str = "UNBOUNDED";  break;
        case SolveStatus::INFEASIBLE: status_str = "INFEASIBLE"; break;
    }
    write("Result: " + status_str);
    {
        std::ostringstream ss;
        ss << "  x          = " << result.x.transpose();
        write(ss.str());
    }
    {
        std::ostringstream ss;
        ss << "  z          = " << result.z;
        write(ss.str());
    }
    {
        std::ostringstream ss;
        ss << "  iterations = " << result.iterations;
        write(ss.str());
    }
}
