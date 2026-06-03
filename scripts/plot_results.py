"""
Generate benchmark visualisations from results/results.csv.
Outputs PNG files to results/.

Usage:
    python scripts/plot_results.py
"""

from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

ROOT = Path(__file__).parent.parent
RESULTS = ROOT / "results"

# ── Palette & style ────────────────────────────────────────────────────────────
LANG_ORDER = ["C++", "Rust", "Julia", "Python", "TypeScript", "Go"]
LANG_COLORS = {
    "C++":        "#2196F3",   # blue
    "Rust":       "#FF5722",   # deep orange
    "Julia":      "#9C27B0",   # purple
    "Python":     "#4CAF50",   # green
    "TypeScript": "#FFC107",   # amber
    "Go":         "#00BCD4",   # cyan
}
PROB_ORDER  = ["small", "medium", "large", "xlarge"]
PROB_LABELS = {"small": "small\n(m=10, n=20)", "medium": "medium\n(m=50, n=100)",
               "large": "large\n(m=100, n=200)", "xlarge": "xlarge\n(m=200, n=400)"}

plt.rcParams.update({
    "font.family":      "sans-serif",
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.color":       "#E0E0E0",
    "grid.linewidth":   0.8,
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  10,
    "figure.dpi":       150,
})


# ── Load data ──────────────────────────────────────────────────────────────────
def load() -> dict[tuple[str, str], dict]:
    data: dict[tuple[str, str], dict] = {}
    with open(RESULTS / "results.csv") as f:
        for row in csv.DictReader(f):
            if row["status"] != "OK":
                continue
            data[(row["language"], row["problem"])] = {
                "median_ms":   float(row["median_ms"]),
                "std_ms":      float(row["std_ms"]),
                "peak_rss_mb": float(row["peak_rss_mb"]),
                "iterations":  int(row["iterations"]),
            }
    return data


# ── Plot 1: Grouped bar — solve time per problem ───────────────────────────────
def plot_solve_times(data: dict) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=False)
    fig.suptitle("Solve Time — Revised Simplex Benchmark", fontsize=14, fontweight="bold", y=1.02)

    for ax, prob in zip(axes, PROB_ORDER):
        langs  = [l for l in LANG_ORDER if (l, prob) in data]
        values = [data[(l, prob)]["median_ms"] for l in langs]
        errs   = [data[(l, prob)]["std_ms"]    for l in langs]
        colors = [LANG_COLORS[l] for l in langs]

        bars = ax.bar(langs, values, yerr=errs, color=colors,
                      capsize=4, width=0.6, error_kw={"linewidth": 1.2})

        # Value labels on bars
        for bar, val in zip(bars, values):
            label = f"{val:.0f}" if val >= 1 else f"{val:.2f}"
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(errs) * 0.05,
                    label, ha="center", va="bottom", fontsize=8, fontweight="bold")

        ax.set_title(PROB_LABELS[prob], fontsize=10)
        ax.set_ylabel("Median time (ms)" if ax == axes[0] else "")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=30)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f"{x:.0f}" if x >= 1 else f"{x:.2f}"))

    fig.tight_layout()
    out = RESULTS / "solve_times.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out.name}")


# ── Plot 2: Memory usage (peak RSS) ───────────────────────────────────────────
def plot_memory(data: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Peak RSS — Memory Usage", fontsize=14, fontweight="bold")

    # Left: bar per language (xlarge problem)
    ax = axes[0]
    prob = "xlarge"
    langs  = [l for l in LANG_ORDER if (l, prob) in data]
    values = [data[(l, prob)]["peak_rss_mb"] for l in langs]
    colors = [LANG_COLORS[l] for l in langs]
    bars = ax.bar(langs, values, color=colors, width=0.6)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                f"{val:.0f} MB", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_title("xlarge problem (m=200, n=400)", fontsize=11)
    ax.set_ylabel("Peak RSS (MB)")
    ax.tick_params(axis="x", rotation=30)

    # Right: grouped lines — RSS vs problem size per language
    ax = axes[1]
    prob_sizes = [10 * 20, 50 * 100, 100 * 200, 200 * 400]   # m*n as proxy for size
    for lang in LANG_ORDER:
        vals = [data.get((lang, p), {}).get("peak_rss_mb") for p in PROB_ORDER]
        if any(v is None for v in vals):
            continue
        ax.plot(PROB_ORDER, vals, marker="o", label=lang,
                color=LANG_COLORS[lang], linewidth=2, markersize=7)
    ax.set_title("RSS across problem sizes", fontsize=11)
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_xlabel("Problem size")
    ax.legend(loc="upper left")

    fig.tight_layout()
    out = RESULTS / "memory.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out.name}")


# ── Plot 3: Log-log scaling ────────────────────────────────────────────────────
def plot_scaling(data: dict) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_title("Solve Time Scaling (log–log)", fontsize=14, fontweight="bold")

    # Use m as problem size proxy
    sizes = {"small": 10, "medium": 50, "large": 100, "xlarge": 200}

    for lang in LANG_ORDER:
        xs, ys = [], []
        for prob in PROB_ORDER:
            entry = data.get((lang, prob))
            if entry:
                xs.append(sizes[prob])
                ys.append(entry["median_ms"])
        if not xs:
            continue
        ax.plot(xs, ys, marker="o", label=lang,
                color=LANG_COLORS[lang], linewidth=2.5, markersize=8)

        # Fit slope (only if 3+ points)
        if len(xs) >= 3:
            log_x = np.log10(xs)
            log_y = np.log10(ys)
            slope = np.polyfit(log_x, log_y, 1)[0]
            ax.annotate(f"  ≈ O(n^{slope:.1f})",
                        xy=(xs[-1], ys[-1]),
                        xytext=(xs[-1] * 1.05, ys[-1]),
                        fontsize=8, color=LANG_COLORS[lang], va="center")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Problem size m (rows)")
    ax.set_ylabel("Median solve time (ms)")
    ax.set_xticks([10, 50, 100, 200])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.legend(loc="upper left", framealpha=0.9)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    fig.tight_layout()
    out = RESULTS / "scaling.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out.name}")


# ── Plot 4: Heatmap — relative slowdown vs fastest ────────────────────────────
def plot_heatmap(data: dict) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title("Relative Slowdown vs Fastest\n(1× = fastest language for that problem)",
                 fontsize=13, fontweight="bold")

    matrix = []
    for lang in LANG_ORDER:
        row = []
        for prob in PROB_ORDER:
            entry = data.get((lang, prob))
            row.append(entry["median_ms"] if entry else float("nan"))
        matrix.append(row)

    matrix = np.array(matrix, dtype=float)
    # Normalise each column by its minimum
    col_min = np.nanmin(matrix, axis=0)
    rel = matrix / col_min

    im = ax.imshow(rel, cmap="RdYlGn_r", aspect="auto", vmin=1, vmax=20)
    plt.colorbar(im, ax=ax, label="× slower than fastest")

    ax.set_xticks(range(len(PROB_ORDER)))
    ax.set_xticklabels([PROB_LABELS[p].replace("\n", " ") for p in PROB_ORDER], fontsize=9)
    ax.set_yticks(range(len(LANG_ORDER)))
    ax.set_yticklabels(LANG_ORDER)

    for i, lang in enumerate(LANG_ORDER):
        for j, prob in enumerate(PROB_ORDER):
            val = rel[i, j]
            txt = f"{val:.1f}×" if not np.isnan(val) else "—"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="white" if val > 8 else "black")

    fig.tight_layout()
    out = RESULTS / "heatmap.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out.name}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main() -> None:
    data = load()
    print("Generating plots...")
    plot_solve_times(data)
    plot_memory(data)
    plot_scaling(data)
    plot_heatmap(data)
    print("Done.")


if __name__ == "__main__":
    main()
