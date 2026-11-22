"""Solve the product-reference KL objective on the grid shared by Figures 4 and 5."""
import os
for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import csv
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.special import logsumexp
from inputs import grids, load_inputs, prepare_inputs

HERE = Path(__file__).resolve().parent
TOL = 1e-10


def control_potential(cost, f, eps, rho, log_a):
    return -(rho / (rho + eps)) * eps * logsumexp(
        log_a[:, None] + (f[:, None] - cost) / eps, axis=0)


def log_sinkhorn(cost, f, eps, rho, log_a, log_b):
    factor = -eps * rho / (rho + eps)
    for n in range(1, 2001):
        g = control_potential(cost, f, eps, rho, log_a)
        new_f = factor * logsumexp(log_b + (g - cost) / eps, axis=1)
        change = np.max(np.abs(new_f - f))
        f = new_f
        if change < 2e-13:
            return f, n
    raise RuntimeError(f"Log Sinkhorn failed at {eps=}, {rho=}")


def newton(cost, f, eps, rho, a, b):
    # With f fixed, optimize the control potential analytically. Each column
    # then has mass s and conditional probabilities p over treated subjects.
    # The remaining convex objective has gradient equal to the row mismatch.
    log_a = np.log(a)

    def evaluate(f, derivatives):
        scores = (f[:, None] - cost) / eps + log_a[:, None]
        peak = scores.max(axis=0)
        p = np.exp(scores - peak)
        sums = p.sum(axis=0)
        s = b * np.exp((eps / (rho + eps)) * (peak + np.log(sums)))
        h = a * np.exp(-f / rho)
        objective = rho * h.sum() + (rho + eps) * s.sum()
        if not derivatives:
            return objective
        p /= sums
        rows = p @ s
        weighted = p * np.sqrt(s)
        hessian = np.diag(rows / eps + h / rho) - (
            rho / (eps * (rho + eps))) * (weighted @ weighted.T)
        return objective, rows - h, hessian, rho * np.max(np.abs(np.log(rows/h)))

    for n in range(1, 121):
        objective, gradient, hessian, residual = evaluate(f, True)
        if residual < 1e-11:
            return f, n
        # Equilibrate before solving. In the small-penalty regime the fitted
        # row masses (and therefore the Hessian diagonals) span many decades.
        scale = np.sqrt(hessian.diagonal())
        system = (hessian / scale[:, None]) / scale
        step = solve(system, -gradient / scale, assume_a="pos", check_finite=False) / scale
        descent = gradient @ step
        rate = 1.0
        for _ in range(45):
            candidate = f + rate * step
            if evaluate(candidate, False) <= objective + 1e-4 * rate * descent + 8e-15 * abs(objective):
                f = candidate
                break
            rate *= 0.5
        else:
            raise RuntimeError(f"Newton line search failed at {eps=}, {rho=}")
    raise RuntimeError(f"Newton failed at {eps=}, {rho=}")


def median_controls_95(log_plan, log_rows):
    # A finite log weight represents positive mass even when exp rounds to zero.
    if not np.all(np.isfinite(log_plan)):
        raise ValueError("Expected a finite log weight for each empirical pair")

    # Rank within each row and accumulate conditional mass, largest weight first.
    conditional = np.exp(log_plan - log_rows[:, None])
    ranked = np.sort(conditional, axis=1)[:, ::-1]
    cumulative = np.cumsum(ranked, axis=1)
    if not np.allclose(cumulative[:, -1], 1, rtol=0, atol=1e-11):
        raise ValueError("Conditional row weights must sum to one")
    controls_95 = 1 + np.sum(cumulative < 0.95, axis=1)
    return np.median(controls_95)


def certify(cost, f, eps, rho, a, b, treated, controls):
    log_a, log_b = np.log(a), np.log(b)
    g = control_potential(cost, f, eps, rho, log_a)
    log_plan = log_a[:, None] + log_b + (f[:, None] + g - cost) / eps
    log_rows = logsumexp(log_plan, axis=1)
    log_cols = logsumexp(log_plan, axis=0)
    # The cellwise gradient of the primal KL objective equals u[i] + v[j].
    u = f + rho * (log_rows - log_a)
    v = g + rho * (log_cols - log_b)
    residual = max(abs(u.min() + v.min()), abs(u.max() + v.max()))
    if not np.isfinite(residual) or residual > TOL:
        raise RuntimeError(f"First-order residual {residual} at {eps=}, {rho=}")
    log_mass = logsumexp(log_rows)
    treated_mean = np.exp(log_rows - log_mass) @ treated
    control_mean = np.exp(log_cols - log_mass) @ controls
    return dict(estimate=treated_mean-control_mean, treated_mean=treated_mean,
                control_mean=control_mean, transport_mass=np.exp(log_mass), max_abs_foc=residual,
                median_controls_95=median_controls_95(log_plan, log_rows))


def solve_row(task):
    rho, eps_values = task
    cost, a, b, yt, yc, benchmark = load_inputs()
    log_a, log_b = np.log(a), np.log(b)
    f = np.zeros_like(a)
    records = []
    for eps in sorted(eps_values, reverse=True):
        if rho <= 2 * eps:
            f, iterations = log_sinkhorn(cost, f, eps, rho, log_a, log_b)
            method = "log_sinkhorn"
        else:
            f, iterations = newton(cost, f, eps, rho, a, b)
            method = "semi_dual_newton"
        stats = certify(cost, f, eps, rho, a, b, yt, yc)
        records.append(dict(epsilon=eps, rho=rho, **stats,
                            experimental_benchmark=benchmark, method=method, iterations=iterations))
    return records


def run(jobs):
    prepare_inputs()
    epsilons, rhos = grids()
    tasks = [(rho, epsilons) for rho in rhos]
    rows = []
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for records in pool.map(solve_row, tasks):
            rows.extend(records)
            print(f"Certified {len(rows)}/{len(epsilons) * len(rhos)} estimates", flush=True)
    rows.sort(key=lambda row: (row["epsilon"], row["rho"]))
    selected = min(rows,
                   key=lambda r: abs(r["estimate"] - r["experimental_benchmark"]))
    for row in rows:
        row["is_selected_table1"] = int(row is selected)
    output = HERE / "results" / "converged_grid.csv"
    output.parent.mkdir(exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=output.parent, newline="", delete=False) as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    Path(stream.name).replace(output)
    print(f"Selected UOT contrast: {selected['estimate']:.2f} at "
          f"epsilon={selected['epsilon']}, rho={selected['rho']}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=2)
    run(parser.parse_args().jobs)
