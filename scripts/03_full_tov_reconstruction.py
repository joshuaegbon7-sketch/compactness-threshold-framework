import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

G = 6.67430e-11
c = 2.99792458e8
M_sun = 1.98847e30
m_n = 1.6749275e-27
MEV_TO_J = 1.602176634e-13
FM3_TO_M3 = 1.0e45
BASE_DIR = "data"
OUT_DIR = "outputs/validation"
os.makedirs(OUT_DIR, exist_ok=True)
EOS_LIST = ["APR", "SLY4", "DD2"]

def load_numeric_file(path):
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric rows found in {path}")
    lengths = [len(row) for row in rows]
    common_length = max(set(lengths), key=lengths.count)
    rows = [row for row in rows if len(row) == common_length]
    return np.array(rows, dtype=float)

def load_official_mr(eos_name):
    data = load_numeric_file(os.path.join(BASE_DIR, eos_name, "eos.mr"))
    R_km, M_msun = data[:, 0], data[:, 1]
    idx = np.argmax(M_msun)
    return R_km, M_msun, {"Mmax": M_msun[idx], "Rmax": R_km[idx]}

def load_compose_eos(eos_name):
    folder = os.path.join(BASE_DIR, eos_name)
    with open(os.path.join(folder, "eos.nb"), "r", errors="ignore") as f:
        nb_lines = [line.strip() for line in f if line.strip()]
    nb = np.array([float(x) for x in nb_lines[2:]], dtype=float)
    thermo_rows = []
    with open(os.path.join(folder, "eos.thermo"), "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                thermo_rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    if not thermo_rows:
        raise ValueError(f"{eos_name}: no numeric rows found in eos.thermo")
    lengths = [len(row) for row in thermo_rows]
    common_length = max(set(lengths), key=lengths.count)
    thermo = np.array([row for row in thermo_rows if len(row) == common_length], dtype=float)
    if len(thermo) == len(nb) + 1:
        thermo = thermo[1:]
    if len(nb) != len(thermo):
        raise ValueError(f"{eos_name}: length mismatch nb={len(nb)}, thermo={len(thermo)}")
    Q_pressure = thermo[:, 3]
    Q_energy = thermo[:, 9]
    nb_m3 = nb * FM3_TO_M3
    P = Q_pressure * nb_m3 * MEV_TO_J
    eps = (1.0 + Q_energy) * nb_m3 * m_n * c**2
    mask = np.isfinite(P) & np.isfinite(eps) & (P > 0) & (eps > 0)
    return P[mask], eps[mask], nb[mask]

def build_eps_of_p(P, eps):
    mask = np.isfinite(P) & np.isfinite(eps) & (P > 0) & (eps > 0)
    P, eps = P[mask], eps[mask]
    order = np.argsort(P)
    P, eps = P[order], eps[order]
    P, idx = np.unique(P, return_index=True)
    eps = eps[idx]
    if len(P) < 20:
        return None, None, None
    interp = interp1d(np.log(P), np.log(eps), bounds_error=False, fill_value="extrapolate")
    def eps_of_p(Pval):
        Pval = np.maximum(Pval, P[0])
        return np.exp(interp(np.log(Pval)))
    return eps_of_p, P[0], P[-1]

def tov_rhs(r, y, eps_of_p):
    m, P = y
    if P <= 0.0 or r <= 0.0:
        return [0.0, 0.0]
    eps = eps_of_p(P)
    compact_factor = 1.0 - 2.0 * G * m / (r * c**2)
    if compact_factor <= 0.0:
        return [0.0, 0.0]
    dm_dr = 4.0 * np.pi * r**2 * eps / c**2
    dP_dr = -G * (eps + P) / c**2
    dP_dr *= m + 4.0 * np.pi * r**3 * P / c**2
    dP_dr /= r**2 * compact_factor
    return [dm_dr, dP_dr]

def integrate_star(Pc, eps_of_p, P_surface):
    r0 = 1.0
    eps_c = eps_of_p(Pc)
    m0 = 4.0 * np.pi * r0**3 * eps_c / (3.0 * c**2)
    def surface_event(r, y):
        return y[1] - P_surface
    surface_event.terminal = True
    surface_event.direction = -1
    sol = solve_ivp(lambda r, y: tov_rhs(r, y, eps_of_p), (r0, 2.0e5), [m0, Pc], events=surface_event, max_step=100.0, rtol=1e-5, atol=1e-8)
    if len(sol.t_events[0]) == 0:
        return None
    return {"R_km": sol.t_events[0][0] / 1000.0, "M_Msun": sol.y_events[0][0][0] / M_sun}

def run_tov_sequence(P, eps, n_central=100):
    eps_of_p, Pmin, Pmax = build_eps_of_p(P, eps)
    if eps_of_p is None:
        return None
    central_pressures = np.logspace(np.log10(Pmin * 100.0), np.log10(Pmax * 0.95), n_central)
    rows = []
    for Pc in central_pressures:
        result = integrate_star(Pc, eps_of_p, Pmin * 1.001)
        if result is not None and 0.5 < result["R_km"] < 120 and 0.001 < result["M_Msun"] < 5:
            rows.append(result)
    return pd.DataFrame(rows) if len(rows) >= 5 else None

summary_rows = []
for eos_name in EOS_LIST:
    print("=" * 70)
    print(f"Running thermo--TOV reconstruction for {eos_name}")
    _, _, official_summary = load_official_mr(eos_name)
    P, eps, _ = load_compose_eos(eos_name)
    seq = run_tov_sequence(P, eps)
    if seq is None:
        print(f"{eos_name}: reconstruction failed.")
        continue
    best = seq.loc[seq["M_Msun"].idxmax()]
    row = {
        "EOS": eos_name,
        "Mmax_reconstructed_Msun": best["M_Msun"],
        "Rmax_reconstructed_km": best["R_km"],
        "Mmax_official_Msun": official_summary["Mmax"],
        "Rmax_official_km": official_summary["Rmax"],
    }
    row["M_frac_error"] = abs(row["Mmax_reconstructed_Msun"] - row["Mmax_official_Msun"]) / row["Mmax_official_Msun"]
    row["R_frac_error"] = abs(row["Rmax_reconstructed_km"] - row["Rmax_official_km"]) / row["Rmax_official_km"]
    summary_rows.append(row)
    seq_path = os.path.join(OUT_DIR, f"{eos_name.lower()}_full_tov_reconstructed_sequence.csv")
    seq.to_csv(seq_path, index=False)
    print(f"Saved: {seq_path}")
summary = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUT_DIR, "full_tov_validation_summary_all.csv")
summary.to_csv(summary_path, index=False)
print(summary)
print(f"Saved: {summary_path}")
