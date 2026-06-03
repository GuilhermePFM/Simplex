// Two-phase Revised Simplex with hand-rolled LU factorization (no Eigen).
// Usage: ./bench_mps_lu <data.json> [N=10]

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// Minimal JSON parser (copied verbatim from bench_mps.cpp — no Eigen)
// ---------------------------------------------------------------------------

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::ostringstream ss; ss << f.rdbuf(); return ss.str();
}
static int parse_int_key(const std::string& json, const std::string& key) {
    std::string pat = "\"" + key + "\"";
    std::size_t pos = json.find(pat);
    if (pos == std::string::npos) throw std::runtime_error("Key not found: " + key);
    pos = json.find(':', pos) + 1;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) ++pos;
    return std::stoi(json.substr(pos));
}
static std::vector<double> parse_number_array(const std::string& json,
                                               std::size_t start, std::size_t& end) {
    std::vector<double> vals;
    if (json[start] != '[') throw std::runtime_error("Expected '['");
    std::size_t pos = start + 1;
    while (pos < json.size()) {
        char c = json[pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == ',') { ++pos; continue; }
        if (c == ']') { end = pos + 1; return vals; }
        std::size_t num_end; double v = std::stod(json.substr(pos), &num_end);
        vals.push_back(v); pos += num_end;
    }
    throw std::runtime_error("Unterminated array");
}
static std::vector<std::vector<double>> parse_2d_array(const std::string& json,
                                                        std::size_t start, std::size_t& end) {
    std::vector<std::vector<double>> rows;
    if (json[start] != '[') throw std::runtime_error("Expected '['");
    std::size_t pos = start + 1;
    while (pos < json.size()) {
        char c = json[pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == ',') { ++pos; continue; }
        if (c == ']') { end = pos + 1; return rows; }
        if (c == '[') {
            std::size_t re; rows.push_back(parse_number_array(json, pos, re)); pos = re; continue;
        }
        throw std::runtime_error("Unexpected char in 2D array");
    }
    throw std::runtime_error("Unterminated 2D array");
}
static std::size_t find_value_start(const std::string& json, const std::string& key) {
    std::string pat = "\"" + key + "\"";
    std::size_t pos = json.find(pat);
    if (pos == std::string::npos) throw std::runtime_error("Key not found: " + key);
    pos = json.find(':', pos) + 1;
    while (pos < json.size() && (json[pos]==' '||json[pos]=='\t'||json[pos]=='\n'||json[pos]=='\r'))
        ++pos;
    return pos;
}

struct LpData { int m, n; std::vector<double> A, b, c; };
static LpData load_lp(const std::string& path) {
    std::string json = read_file(path);
    LpData lp; lp.m = parse_int_key(json,"m"); lp.n = parse_int_key(json,"n");
    std::size_t a_start = find_value_start(json,"A"), a_end;
    auto rows = parse_2d_array(json, a_start, a_end);
    lp.A.resize(lp.m * lp.n);
    for (int i = 0; i < lp.m; ++i)
        for (int j = 0; j < lp.n; ++j) lp.A[i*lp.n+j] = rows[i][j];
    std::size_t b_start = find_value_start(json,"b"), b_end;
    lp.b = parse_number_array(json, b_start, b_end);
    std::size_t c_start = find_value_start(json,"c"), c_end;
    lp.c = parse_number_array(json, c_start, c_end);
    return lp;
}

// ---------------------------------------------------------------------------
// LU factorization (row-major flat vector)
// ---------------------------------------------------------------------------

struct LU {
    std::vector<double> m;
    std::vector<int> perm;
    int sz;
};

static LU lu_factor(const std::vector<double>& B, int sz) {
    LU lu; lu.sz = sz; lu.m = B;
    lu.perm.resize(sz); std::iota(lu.perm.begin(), lu.perm.end(), 0);
    for (int k = 0; k < sz; ++k) {
        int r = k; double maxv = std::abs(lu.m[k*sz+k]);
        for (int i = k+1; i < sz; ++i) {
            double v = std::abs(lu.m[i*sz+k]);
            if (v > maxv) { maxv = v; r = i; }
        }
        if (r != k) {
            for (int j = 0; j < sz; ++j) std::swap(lu.m[k*sz+j], lu.m[r*sz+j]);
            std::swap(lu.perm[k], lu.perm[r]);
        }
        if (std::abs(lu.m[k*sz+k]) < 1e-14) continue;
        for (int i = k+1; i < sz; ++i) {
            lu.m[i*sz+k] /= lu.m[k*sz+k];
            for (int j = k+1; j < sz; ++j)
                lu.m[i*sz+j] -= lu.m[i*sz+k] * lu.m[k*sz+j];
        }
    }
    return lu;
}

static std::vector<double> lu_solve_vec(const LU& lu, const std::vector<double>& b) {
    int sz = lu.sz;
    std::vector<double> y(sz);
    for (int i = 0; i < sz; ++i) y[i] = b[lu.perm[i]];
    for (int i = 1; i < sz; ++i)
        for (int k = 0; k < i; ++k) y[i] -= lu.m[i*sz+k]*y[k];
    std::vector<double> x = y;
    for (int i = sz-1; i >= 0; --i) {
        for (int k = i+1; k < sz; ++k) x[i] -= lu.m[i*sz+k]*x[k];
        if (std::abs(lu.m[i*sz+i]) > 1e-14) x[i] /= lu.m[i*sz+i];
    }
    return x;
}

// Solve B·X = N; N is m×nb row-major; returns m×nb row-major.
static std::vector<double> lu_solve_mat(const LU& lu, const std::vector<double>& N, int nb) {
    int m = lu.sz;
    std::vector<double> res(m*nb), col(m);
    for (int j = 0; j < nb; ++j) {
        for (int i = 0; i < m; ++i) col[i] = N[i*nb+j];
        auto x = lu_solve_vec(lu, col);
        for (int i = 0; i < m; ++i) res[i*nb+j] = x[i];
    }
    return res;
}

// ---------------------------------------------------------------------------
// Matrix helpers
// ---------------------------------------------------------------------------

static std::vector<double> cols_rm(const std::vector<double>& A, int m, int n,
                                    const std::vector<int>& idx) {
    int k = (int)idx.size();
    std::vector<double> out(m*k);
    for (int jj = 0; jj < k; ++jj) {
        int col = idx[jj];
        for (int i = 0; i < m; ++i) out[i*k+jj] = A[i*n+col];
    }
    return out;
}

static double vdot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0; for (int i = 0; i < (int)a.size(); ++i) s += a[i]*b[i]; return s;
}

static std::vector<double> gather(const std::vector<double>& a, const std::vector<int>& idx) {
    std::vector<double> out(idx.size());
    for (int i = 0; i < (int)idx.size(); ++i) out[i] = a[idx[i]];
    return out;
}

// cr = cN - BinvN'·cB; BinvN is m×nb row-major
static std::vector<double> reduced_costs(const std::vector<double>& cN,
                                          const std::vector<double>& cB,
                                          const std::vector<double>& BinvN, int m, int nb) {
    std::vector<double> cr = cN;
    for (int j = 0; j < nb; ++j)
        for (int i = 0; i < m; ++i) cr[j] -= BinvN[i*nb+j]*cB[i];
    return cr;
}

static std::vector<double> dense_col(const std::vector<double>& BinvN, int m, int nb, int j) {
    std::vector<double> col(m);
    for (int i = 0; i < m; ++i) col[i] = BinvN[i*nb+j];
    return col;
}

// ---------------------------------------------------------------------------
// Pivot selection
// ---------------------------------------------------------------------------

static int entering_index(const std::vector<double>& cr, double tol) {
    int best = -1; double bv = tol;
    for (int j = 0; j < (int)cr.size(); ++j)
        if (cr[j] > bv) { bv = cr[j]; best = j; }
    return best;
}

static int leaving_index(const std::vector<double>& xb, const std::vector<double>& d,
                          const std::vector<double>& BinvN, int m, int nb, double tol) {
    std::vector<int> cands;
    for (int i = 0; i < m; ++i) if (d[i] > tol) cands.push_back(i);
    if (cands.empty()) return -1;
    int best = cands[0];
    for (int ii = 1; ii < (int)cands.size(); ++ii) {
        int i = cands[ii];
        double ri = xb[i]/d[i], rb = xb[best]/d[best];
        if (ri < rb - tol) { best = i; }
        else if (std::abs(ri-rb) <= tol) {
            for (int k = 0; k < nb; ++k) {
                double vi = BinvN[i*nb+k]/d[i], vb = BinvN[best*nb+k]/d[best];
                if (vi < vb-1e-14) { best = i; break; }
                else if (vi > vb+1e-14) break;
            }
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// Phase 1
// ---------------------------------------------------------------------------

static void make_b_nonneg(std::vector<double>& A, std::vector<double>& b, int m, int n) {
    for (int i = 0; i < m; ++i) if (b[i] < 0) {
        b[i] = -b[i];
        for (int j = 0; j < n; ++j) A[i*n+j] = -A[i*n+j];
    }
}

struct BFS { std::vector<int> basic, nonbasic; bool ok; };

static void fix_arts(std::vector<double>& Aw, std::vector<int>& basic,
                     std::vector<int>& nonbasic, int m, int nw, int artStart) {
    for (bool changed = true; changed;) {
        changed = false;
        for (int pos = 0; pos < (int)basic.size(); ++pos) {
            if (basic[pos] < artStart) continue;
            auto B = cols_rm(Aw, m, nw, basic);
            auto lu = lu_factor(B, m);
            auto BinvAw = lu_solve_mat(lu, Aw, nw);
            int swapJ = -1;
            for (int j = 0; j < (int)nonbasic.size(); ++j)
                if (nonbasic[j] < artStart && std::abs(BinvAw[pos*nw+nonbasic[j]]) > 1e-10)
                    { swapJ = j; break; }
            if (swapJ < 0) {
                for (int j = 0; j < nw; ++j) Aw[pos*nw+j] = 0.0;
            } else {
                std::swap(basic[pos], nonbasic[swapJ]);
                changed = true; break;
            }
        }
    }
}

static BFS phase1(const std::vector<double>& A_orig, const std::vector<double>& b_orig,
                  int m, int n, int maxit=10000, double tol=1e-8) {
    std::vector<double> A = A_orig, b = b_orig;
    make_b_nonneg(A, b, m, n);
    int nw = n + m;
    std::vector<double> Aw(m*nw, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) Aw[i*nw+j] = A[i*n+j];
        Aw[i*nw+n+i] = 1.0;
    }
    std::vector<double> cw(nw, 0.0);
    for (int i = n; i < nw; ++i) cw[i] = -1.0;

    std::vector<int> basic(m), nonbasic(n);
    std::iota(basic.begin(), basic.end(), n);
    std::iota(nonbasic.begin(), nonbasic.end(), 0);

    for (int it = 0; it < maxit; ++it) {
        int nb = (int)nonbasic.size();
        auto B = cols_rm(Aw, m, nw, basic);
        auto N = cols_rm(Aw, m, nw, nonbasic);
        auto lu = lu_factor(B, m);
        auto xb = lu_solve_vec(lu, b);
        auto BinvN = lu_solve_mat(lu, N, nb);
        auto cb = gather(cw, basic);
        auto cr = reduced_costs(gather(cw, nonbasic), cb, BinvN, m, nb);
        double z = vdot(cb, xb);

        int j = entering_index(cr, tol);
        if (j < 0) {
            if (z < -tol) return {{},{},false};
            fix_arts(Aw, basic, nonbasic, m, nw, n);
            std::vector<int> bas, nonbas;
            for (int v : basic)   if (v < n) bas.push_back(v);
            for (int v : nonbasic) if (v < n) nonbas.push_back(v);
            return {bas, nonbas, true};
        }
        auto d = dense_col(BinvN, m, nb, j);
        int i = leaving_index(xb, d, BinvN, m, nb, tol);
        if (i < 0) throw std::runtime_error("Phase 1 unexpectedly unbounded");
        std::swap(basic[i], nonbasic[j]);
    }
    throw std::runtime_error("Phase 1 did not converge");
}

// ---------------------------------------------------------------------------
// Phase 2
// ---------------------------------------------------------------------------

static std::pair<double,int> phase2(const std::vector<double>& A, const std::vector<double>& b,
                                     const std::vector<double>& c, int m, int n,
                                     std::vector<int> basic, std::vector<int> nonbasic,
                                     int maxit=10000, double tol=1e-10) {
    for (int it = 1; it <= maxit; ++it) {
        int nb = (int)nonbasic.size();
        auto B = cols_rm(A, m, n, basic);
        auto N = cols_rm(A, m, n, nonbasic);
        auto lu = lu_factor(B, m);
        auto xb = lu_solve_vec(lu, b);
        auto BinvN = lu_solve_mat(lu, N, nb);
        auto cb = gather(c, basic);
        auto cr = reduced_costs(gather(c, nonbasic), cb, BinvN, m, nb);
        double z = vdot(cb, xb);

        int j = entering_index(cr, tol);
        if (j < 0) return {z, it};
        auto d = dense_col(BinvN, m, nb, j);
        int i = leaving_index(xb, d, BinvN, m, nb, tol);
        if (i < 0) return {std::numeric_limits<double>::infinity(), it};
        std::swap(basic[i], nonbasic[j]);
    }
    throw std::runtime_error("Phase 2 did not converge");
}

static std::pair<double,int> simplex_solve(const LpData& lp) {
    int m = lp.m, n = lp.n;
    auto bfs = phase1(lp.A, lp.b, m, n);
    if (!bfs.ok) return {0.0, 0};
    // fill any unaccounted columns into nonbasic
    std::unordered_map<int,bool> inSet;
    for (int v : bfs.basic)    inSet[v] = true;
    for (int v : bfs.nonbasic) inSet[v] = true;
    for (int j = 0; j < n; ++j)
        if (!inSet[j]) bfs.nonbasic.push_back(j);
    return phase2(lp.A, lp.b, lp.c, m, n, bfs.basic, bfs.nonbasic);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cerr << "Usage: bench_mps_lu <data.json> [N]\n"; return 1; }
    int n_runs = (argc >= 3) ? std::atoi(argv[2]) : 10;
    LpData lp = load_lp(argv[1]);

    using Clock = std::chrono::high_resolution_clock;
    using Dur   = std::chrono::duration<double, std::milli>;

    double last_z = 0.0; int last_it = 0;
    std::vector<double> times_ms; times_ms.reserve(n_runs);

    for (int r = 0; r < n_runs; ++r) {
        auto t0 = Clock::now();
        auto [z, it] = simplex_solve(lp);
        times_ms.push_back(Dur(Clock::now() - t0).count());
        last_z = z; last_it = it;
    }

    std::cout << "{\"z\": " << std::fixed; std::cout.precision(4);
    std::cout << last_z << ", \"iterations\": " << last_it << ", \"times_ms\": [";
    std::cout.precision(6);
    for (int i = 0; i < n_runs; ++i) { if (i) std::cout << ", "; std::cout << times_ms[i]; }
    std::cout << "]}\n";
    return 0;
}
