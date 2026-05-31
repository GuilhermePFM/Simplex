@enum LogLevel SILENT INFO DEBUG

struct SimplexLogger
    level  :: LogLevel
    stream :: IO
end
SimplexLogger() = SimplexLogger(SILENT, devnull)
SimplexLogger(level::LogLevel) = SimplexLogger(level, stdout)
SimplexLogger(path::AbstractString, level::LogLevel = DEBUG) = SimplexLogger(level, open(path, "w"))

function log_problem(logger::SimplexLogger, lp::LinearProgram)
    logger.level == SILENT && return
    _write(logger, "=======================")
    _write(logger, " Revised Simplex Solver")
    _write(logger, "=======================")
    _write(logger, "A =\n$(lp.A)")
    _write(logger, "b = $(lp.b)")
    _write(logger, "c = $(lp.c)")
    _write(logger, "")
end

function log_phase(logger::SimplexLogger, phase::Int)
    logger.level == SILENT && return
    header = phase == 1 ? "Phase 1 — Finding initial BFS" :
                          "Phase 2 — Optimizing"
    _write(logger, header)
    _write(logger, repeat("─", length(header)))
end

function log_iteration(logger::SimplexLogger, it::Int, basis::BasisState,
                        xb::AbstractVector{Float64}, z::Float64)
    logger.level != DEBUG && return
    _write(logger, "  iter $it:")
    _write(logger, "    xb       = $xb")
    _write(logger, "    basic    = $(basis.basic)")
    _write(logger, "    nonbasic = $(basis.nonbasic)")
    _write(logger, "    z        = $z")
    _write(logger, "")
end

function log_result(logger::SimplexLogger, result::SimplexResult)
    logger.level == SILENT && return
    _write(logger, "")
    _write(logger, "Result: $(result.status)")
    _write(logger, "  x          = $(result.x)")
    _write(logger, "  z          = $(result.z)")
    _write(logger, "  iterations = $(result.iterations)")
end

function close_log(logger::SimplexLogger)
    logger.stream !== stdout && logger.stream !== devnull && close(logger.stream)
end

_write(l::SimplexLogger, s::AbstractString) = println(l.stream, s)
