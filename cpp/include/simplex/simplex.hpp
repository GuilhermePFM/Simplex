#pragma once

#include "types.hpp"
#include "pivot.hpp"
#include "preprocess.hpp"
#include "logger.hpp"
#include "phase1.hpp"
#include "phase2.hpp"

SimplexResult solve(const LinearProgram& lp,
                    const PivotRule&     rule      = LargestCoefficient{},
                    LogLevel             verbosity = LogLevel::SILENT);

SimplexResult solve_phase2(const LinearProgram& lp,
                            const BasisState&    init_basis,
                            const PivotRule&     rule      = LargestCoefficient{},
                            LogLevel             verbosity = LogLevel::SILENT);
