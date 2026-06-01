#pragma once
#include "types.hpp"
#include "pivot.hpp"
#include "logger.hpp"

SimplexResult
phase2(const LinearProgram& lp, const BasisState& init_basis,
       const PivotRule& rule, const SimplexLogger& logger,
       int maxit = 1000, double tol = 1e-10);
