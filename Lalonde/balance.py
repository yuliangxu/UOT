"""Selected-fit covariate balance on the sample and scale used for the LaLonde cost."""
import numpy as np


# Column indices in the aligned NSW/PSID files (treatment is column zero).
COVARIATES = (
    (1, "age", "Age (years)", "years"),
    (2, "educ", "Education (years)", "years"),
    (7, "re75", "1975 earnings (\\$)", "dollars"),
    (3, "black", "Black (\\%)", "proportion"),
    (4, "hispanic", "Hispanic (\\%)", "proportion"),
    (5, "married", "Married (\\%)", "proportion"),
    (6, "nodegree", "No high-school diploma (\\%)", "proportion"),
)


def format_mean(value, unit):
    if unit == "proportion":
        percentage = 100 * value
        return "$<0.1$" if 0 < percentage < 0.1 else f"{percentage:.1f}"
    return f"{value:,.2f}"


def selected_diagnostics(treated, controls, weights, selected):
    treated_weights, control_weights = weights
    setting = {"epsilon": selected["epsilon"], "rho": selected["rho"]}
    concentration = []
    for group, w in (("NSW treated", treated_weights), ("PSID control rows", control_weights)):
        if not np.all(np.isfinite(w)) or np.any(w < 0) or not np.isclose(w.sum(), 1, atol=1e-10, rtol=0):
            raise ValueError(f"Invalid marginal weights for {group}")
        concentration.append(dict(**setting, group=group, observations=len(w),
                                  effective_sample_size=1 / (w @ w), largest_weight=w.max(),
                                  top_10_weight=np.sort(w)[-10:].sum()))

    balance = []
    tex = [
        r"\begin{tabular}{@{}lrrrrrr@{}}",
        r"\toprule",
        r"& \multicolumn{3}{c}{Before weighting} & \multicolumn{3}{c}{After UOT weighting}\\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
        r"Covariate & Treated & Controls & SMD & Treated & Controls & SMD\\",
        r"\midrule",
    ]
    for column, name, label, unit in COVARIATES:
        x, y = treated[:, column], controls[:, column]
        # Fixed denominator for both SMDs: the sample SD across all original rows.
        pooled_sd = np.concatenate((x, y)).std(ddof=1)
        before_t, before_c = x.mean(), y.mean()
        after_t, after_c = treated_weights @ x, control_weights @ y
        before_smd = (before_t - before_c) / pooled_sd
        after_smd = (after_t - after_c) / pooled_sd
        balance.append(dict(**setting, covariate=name, unit=unit, pooled_sd=pooled_sd,
                            treated_before=before_t, control_before=before_c, smd_before=before_smd,
                            treated_weighted=after_t, control_weighted=after_c, smd_weighted=after_smd))
        # Converting a rounded zero to 0.0 suppresses a misleading minus sign.
        columns = (label, format_mean(before_t, unit), format_mean(before_c, unit),
                   f"{before_smd:.3f}", format_mean(after_t, unit), format_mean(after_c, unit),
                   f"{round(after_smd, 6) or 0.0:.6f}")
        tex.append(" & ".join(columns) + r" \\")
    tex.extend((r"\bottomrule", r"\end{tabular}"))
    return balance, concentration, "\n".join(tex) + "\n"
