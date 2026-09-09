"""
============================================================================
Hausman Test Monte Carlo -- INDEPENDENT PYTHON REIMPLEMENTATION

Independently reproduces the manuscript's core simulation: the full 192,000-
run design (192 scenarios x 1,000 replications). Shares no code with the R
implementation (hausman_full_reproducibility.R); built from the manuscript's
equations directly, using pandas and linearmodels rather than plm.

Feasibility: the aggregate rank condition N x (T-1) >= k holds by at least
two orders of magnitude margin for every one of the 192 architecture x
complexity combinations, so no scenario is excluded a priori.

Output: raw_results_python.csv, one row per replication, matching the
column structure of the R implementation's r_raw_results.csv, so results
from both languages can be pooled or compared directly.
============================================================================
"""
import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS, RandomEffects
from scipy import stats
import time

OUTPUT_DIR = "./results"
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

SEED = 1978
N_REPLICATIONS = 1000  # matches the R implementation; reduce for a quick test run

PANEL_ARCHITECTURES = {
    "Wide Panel": (400, 4),
    "Square": (200, 8),
    "Long Panel": (100, 16),
}
COMPLEXITY_LEVELS = {
    "Simple": 1, "Standard": 3, "Complex": 6, "High-Dimensional": 10,
}
MECHANISMS = ["Random", "Early Exit", "Late Missing", "Cyclical"]
DROPOUT_RATES = [0.10, 0.20, 0.30, 0.40]


def generate_panel(N, T, k, rng):
    ids = np.repeat(np.arange(N), T)
    times = np.tile(np.arange(1, T + 1), N)  # 1-indexed, matches manuscript's t = 1,...,T
    alpha = rng.normal(0, 2, N)              # sd=2 -> Var=4
    alpha_long = np.repeat(alpha, T)
    X = rng.normal(0, 1, size=(N * T, k))
    eps = rng.normal(0, 1.5, N * T)          # sd=1.5 -> Var=2.25
    beta = np.ones(k)
    Y = X @ beta + alpha_long + eps
    df = pd.DataFrame(X, columns=[f"x{i+1}" for i in range(k)])
    df["id"] = ids
    df["time"] = times
    df["y"] = Y
    return df


def apply_missingness(df, mechanism, delta, x_cols, rng):
    d = df.copy()
    T_total = d["time"].nunique()
    N = d["id"].nunique()

    if mechanism == "Random":
        mask = rng.random(len(d)) < delta
        d.loc[mask, ["y"] + x_cols] = np.nan

    elif mechanism == "Early Exit":
        n_drop = max(1, round(2 * N * delta))
        if delta > 0 and n_drop > 0:
            unit_ids = d["id"].unique()
            drop_ids = rng.choice(unit_ids, size=n_drop, replace=False)
            for uid in drop_ids:
                t_exit = rng.integers(2, T_total + 1)
                mask = (d["id"] == uid) & (d["time"] >= t_exit)
                d.loc[mask, ["y"] + x_cols] = np.nan

    elif mechanism == "Late Missing":
        c_val = 2 * delta / (T_total + 1)
        prob = np.minimum(c_val * d["time"], 0.95)
        mask = rng.random(len(d)) < prob
        d.loc[mask, ["y"] + x_cols] = np.nan

    elif mechanism == "Cyclical":
        lo, hi = int(np.floor(T_total / 2)), int(np.ceil(T_total / 2 + 1))
        high_periods = set(range(lo, hi + 1))
        T_high = len(high_periods)
        T_low = T_total - T_high
        p_low = T_total * delta / (3 * T_high + T_low)
        p_high = min(3 * p_low, 0.95)
        probs = d["time"].apply(lambda t: p_high if t in high_periods else p_low)
        mask = rng.random(len(d)) < probs
        d.loc[mask, ["y"] + x_cols] = np.nan

    return d


def hausman_test(df, x_cols):
    d = df.dropna(subset=["y"] + x_cols).set_index(["id", "time"])
    if d.empty or d.reset_index()["id"].nunique() < 2:
        return None, "insufficient_data"
    try:
        exog = d[x_cols]
        fe_res = PanelOLS(d["y"], exog, entity_effects=True).fit()
        re_res = RandomEffects(d["y"], exog).fit()
    except Exception:
        return None, "estimation_error"
    diff = re_res.params - fe_res.params
    var_diff = fe_res.cov - re_res.cov
    try:
        inv_var_diff = np.linalg.inv(var_diff.values)
    except np.linalg.LinAlgError:
        return None, "singular_matrix"
    H = diff.values @ inv_var_diff @ diff.values
    if H < 0:
        # Finite-sample noise can push the estimated variance difference
        # slightly off positive-semi-definite; clip to the boundary (H=0)
        # rather than treat as a computational failure, consistent with
        # standard practice for this known edge case.
        H = 0.0
    p_value = 1 - stats.chi2.cdf(H, len(x_cols))
    return p_value, "success"


def run_simulation():
    # --- Checkpoint / resume setup ---
    # Output path is resolved to an absolute path up front, and a checkpoint
    # file tracks which scenarios are already done. If this script is
    # interrupted (power loss, USB/network drive disconnect, crash, etc.),
    # rerunning it picks up where it left off instead of starting over.
    output_dir_abs = os.path.abspath(OUTPUT_DIR)
    os.makedirs(output_dir_abs, exist_ok=True)
    results_path = os.path.join(output_dir_abs, "raw_results_python.csv")
    checkpoint_path = os.path.join(output_dir_abs, "_checkpoint.txt")

    print(f"Writing results to: {results_path}")
    if os.path.abspath(os.getcwd()).split(os.sep)[0] != output_dir_abs.split(os.sep)[0]:
        print("NOTE: output drive differs from current drive; this is fine, "
              "just confirming the absolute path above is where your data will land.")

    completed_scenarios = set()
    if os.path.exists(checkpoint_path):
        with open(checkpoint_path, "r") as f:
            completed_scenarios = set(line.strip() for line in f if line.strip())
        print(f"Resuming: {len(completed_scenarios)} scenario(s) already completed "
              f"per checkpoint file.")

    file_exists = os.path.exists(results_path) and len(completed_scenarios) > 0

    rng = np.random.default_rng(SEED)
    total_scenarios = len(PANEL_ARCHITECTURES) * len(COMPLEXITY_LEVELS) * len(MECHANISMS) * len(DROPOUT_RATES)
    scenario_num = 0
    start = time.time()

    for panel_name, (N, T) in PANEL_ARCHITECTURES.items():
        for complexity_name, k in COMPLEXITY_LEVELS.items():
            x_cols = [f"x{i+1}" for i in range(k)]
            for mechanism in MECHANISMS:
                for delta in DROPOUT_RATES:
                    scenario_num += 1
                    scenario_key = f"{panel_name}|{complexity_name}|{mechanism}|{delta}"

                    if scenario_key in completed_scenarios:
                        continue  # already done in a prior run; skip but keep RNG state moving on is not needed since skipped scenarios don't consume draws here

                    rows = []
                    for rep in range(N_REPLICATIONS):
                        df = generate_panel(N, T, k, rng)
                        d = apply_missingness(df, mechanism, delta, x_cols, rng)
                        p_value, failure_reason = hausman_test(d, x_cols)
                        specificity = int(p_value > 0.05) if p_value is not None else np.nan
                        rows.append({
                            "p_value": p_value,
                            "failure_reason": failure_reason,
                            "panel_name": panel_name,
                            "complexity_name": complexity_name,
                            "mechanism": mechanism,
                            "delta": delta,
                            "N": N, "T": T, "k": k,
                            "specificity": specificity,
                        })

                    # Append this scenario's results to disk immediately.
                    scenario_df = pd.DataFrame(rows)
                    scenario_df.to_csv(results_path, mode="a", header=not file_exists, index=False)
                    file_exists = True

                    # Record this scenario as done, immediately, on disk.
                    with open(checkpoint_path, "a") as f:
                        f.write(scenario_key + "\n")

                    elapsed = time.time() - start
                    print(f"Scenario {scenario_num}/{total_scenarios} done "
                          f"({elapsed/60:.1f} min elapsed) -- saved to disk")

    print(f"\nDone. All {total_scenarios} scenarios written to {results_path}")
    print("You can delete the _checkpoint.txt file now that the run is complete.")
    return pd.read_csv(results_path)


if __name__ == "__main__":
    cwd = os.path.abspath(os.getcwd())
    if len(cwd) >= 2 and cwd[1] == ':' and cwd[0].upper() not in ('C',):
        print(f"WARNING: running from drive {cwd[0]}:, which may be removable media")
        print("(a USB or external drive). If this drive is disconnected while")
        print("the script is running, the run WILL be interrupted, though progress")
        print("is saved after every scenario and can be resumed by rerunning.")
        print("Recommended: copy this script to your local hard drive (e.g. C:)")
        print("and run it from there instead.\n")
    run_simulation()
