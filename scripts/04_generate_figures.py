import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

G = 6.67430e-11
c = 2.99792458e8
M_sun = 1.98847e30
BASE_DIR = "data"
FIG_DIR = "outputs/figures"
TABLE_DIR = "outputs/tables"
VAL_DIR = "outputs/validation"
os.makedirs(FIG_DIR, exist_ok=True)
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
        raise ValueError(f"No numeric rows in {path}")
    lengths = [len(r) for r in rows]
    common = max(set(lengths), key=lengths.count)
    return np.array([r for r in rows if len(r) == common], dtype=float)

def load_mr(eos):
    data = load_numeric_file(os.path.join(BASE_DIR, eos, "eos.mr"))
    return data[:, 0], data[:, 1]

def compactness(M_msun, R_km):
    return G * (M_msun * M_sun) / ((R_km * 1000.0) * c**2)

def nearest_by_mass(R, M, target=1.4):
    idx = np.argmin(np.abs(M - target))
    return R[idx], M[idx], compactness(M[idx], R[idx])

# 1. Thermo--TOV validation
plt.figure(figsize=(8, 6))
for eos in EOS_LIST:
    R, M = load_mr(eos)
    seq = pd.read_csv(os.path.join(VAL_DIR, f"{eos.lower()}_full_tov_reconstructed_sequence.csv"))
    plt.plot(R, M, linewidth=2.8, label=f"{eos} official")
    plt.plot(seq["R_km"], seq["M_Msun"], linestyle="--", linewidth=2.4, label=f"{eos} reconstructed")
plt.xlabel("Radius R (km)")
plt.ylabel("Mass M ($M_\\odot$)")
plt.title("Thermo--TOV Reconstruction Validation")
plt.xlim(8, 16)
plt.ylim(0.5, 2.7)
plt.grid(True)
plt.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "full_tov_reconstruction_validation_all.png"), dpi=400, bbox_inches="tight")
plt.close()

# 2. Mass-radius vs scaling thresholds
thresholds = pd.read_csv(os.path.join(TABLE_DIR, "scaling_threshold_table.csv"))
plt.figure(figsize=(8, 6))
for eos in EOS_LIST:
    R, M = load_mr(eos)
    plt.plot(R, M, linewidth=2.8, label=eos)
for lam in [0.18, 0.20, 0.22]:
    row = thresholds.iloc[np.argmin(np.abs(thresholds["lambda"] - lam))]
    plt.axhline(row["Mcrit_Msun"], linestyle="--", linewidth=2.2, label=f"Scaling λ={lam:.2f}")
plt.xlabel("Radius R (km)")
plt.ylabel("Mass M ($M_\\odot$)")
plt.title("CompOSE Mass-Radius Curves vs Scaling Thresholds")
plt.xlim(8, 16)
plt.ylim(1.0, 2.6)
plt.grid(True)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "compose_mr_scaling_thresholds_all.png"), dpi=400, bbox_inches="tight")
plt.close()

# 3. Compactness curves
plt.figure(figsize=(8, 6))
for eos in EOS_LIST:
    R, M = load_mr(eos)
    C = compactness(M, R)
    mask = M >= 1.0
    plt.plot(M[mask], C[mask], linewidth=2.8, label=eos)
plt.xlabel("Mass M ($M_\\odot$)")
plt.ylabel("Compactness GM/(Rc$^2$)")
plt.title("Compactness Along CompOSE Mass-Radius Curves")
plt.xlim(1.0, 2.5)
plt.ylim(0.10, 0.34)
plt.grid(True)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "compose_compactness_curves_all.png"), dpi=400, bbox_inches="tight")
plt.close()

# 4. Approx Lambda full range
k2 = 0.08
plt.figure(figsize=(8, 6))
for eos in EOS_LIST:
    R, M = load_mr(eos)
    C = compactness(M, R)
    Lambda = (2.0 / 3.0) * k2 * C**(-5)
    mask = np.isfinite(Lambda) & (M > 0.05)
    plt.semilogy(M[mask], Lambda[mask], linewidth=2.8, label=eos)
plt.xlabel("Mass M ($M_\\odot$)")
plt.ylabel("Approximate $\\Lambda$")
plt.title("Approximate Tidal Deformability")
plt.xlim(0.0, 2.55)
plt.ylim(10, 1e12)
plt.grid(True, which="both")
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "approx_lambda_vs_mass.png"), dpi=400, bbox_inches="tight")
plt.close()

# 5. Approx Lambda NS range
plt.figure(figsize=(8, 6))
for eos in EOS_LIST:
    R, M = load_mr(eos)
    C = compactness(M, R)
    Lambda = (2.0 / 3.0) * k2 * C**(-5)
    mask = (M >= 1.0) & (M <= 2.2) & np.isfinite(Lambda)
    plt.semilogy(M[mask], Lambda[mask], linewidth=2.8, label=eos)
plt.xlabel("Mass M ($M_\\odot$)")
plt.ylabel("Approximate $\\Lambda$")
plt.title("Approximate Tidal Deformability in Neutron-Star Mass Range")
plt.xlim(1.0, 2.25)
plt.ylim(10, 3000)
plt.grid(True, which="both")
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "approx_lambda_ns_range.png"), dpi=400, bbox_inches="tight")
plt.close()

rows = []
for eos in EOS_LIST:
    R, M = load_mr(eos)
    R14, M14, C14 = nearest_by_mass(R, M, target=1.4)
    Lambda14 = (2.0 / 3.0) * k2 * C14**(-5)
    rows.append({"EOS": eos, "M_Msun": M14, "R_km": R14, "Compactness": C14, "Lambda_k2_0p08": Lambda14})
pd.DataFrame(rows).to_csv(os.path.join(TABLE_DIR, "approx_tidal_deformability_summary.csv"), index=False)
print("Generated figures in:", FIG_DIR)
