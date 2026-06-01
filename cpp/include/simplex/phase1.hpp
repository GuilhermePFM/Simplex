#pragma once
#include <tuple>
#include "types.hpp"
#include "pivot.hpp"
#include "logger.hpp"

std::tuple<BasisState, SolveStatus>
phase1(const LinearProgram& lp, const PivotRule& rule, const SimplexLogger& logger,
       int maxit = 1000, double tol = 1e-8);
