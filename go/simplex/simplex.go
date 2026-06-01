package simplex

// Solve runs the two-phase Revised Simplex Method on lp.
// Phase 1 finds an initial BFS; Phase 2 optimizes the objective.
func Solve(lp LinearProgram, opts ...Option) SimplexResult {
	cfg := defaultConfig()
	for _, o := range opts {
		o(&cfg)
	}

	logger := cfg.logger
	logger.LogProblem(lp)

	basis, p1Status := Phase1(lp, cfg.rule, logger, cfg.maxit, cfg.tol1)
	if p1Status == Infeasible {
		logger.Close()
		return SimplexResult{
			X:          make([]float64, lp.N),
			Z:          0.0,
			Status:     Infeasible,
			Iterations: 0,
		}
	}

	result := Phase2(lp, basis, cfg.rule, logger, cfg.maxit, cfg.tol2)
	logger.LogResult(result)
	logger.Close()
	return result
}

// SolvePhase2 runs Phase 2 only, starting from the given initial basis.
// The caller must provide a valid basic feasible solution.
func SolvePhase2(lp LinearProgram, initBasis BasisState, opts ...Option) SimplexResult {
	cfg := defaultConfig()
	for _, o := range opts {
		o(&cfg)
	}

	logger := cfg.logger
	result := Phase2(lp, initBasis, cfg.rule, logger, cfg.maxit, cfg.tol2)
	logger.LogResult(result)
	return result
}

// solveConfig holds configuration for a solve call.
type solveConfig struct {
	rule   PivotRule
	logger SimplexLogger
	maxit  int
	tol1   float64 // Phase 1 tolerance
	tol2   float64 // Phase 2 tolerance
}

func defaultConfig() solveConfig {
	return solveConfig{
		rule:   LargestCoefficient{},
		logger: NewSilentLogger(),
		maxit:  1000,
		tol1:   1e-8,
		tol2:   1e-10,
	}
}

// Option is a functional option for Solve / SolvePhase2.
type Option func(*solveConfig)

// WithRule sets the pivot rule.
func WithRule(r PivotRule) Option {
	return func(c *solveConfig) { c.rule = r }
}

// WithLogger sets the logger.
func WithLogger(l SimplexLogger) Option {
	return func(c *solveConfig) { c.logger = l }
}

// WithMaxIt sets the maximum number of iterations per phase.
func WithMaxIt(n int) Option {
	return func(c *solveConfig) { c.maxit = n }
}
