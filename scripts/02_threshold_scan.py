import os
import numpy as np
import pandas as pd

G = 6.67430e-11
c = 2.99792458e8
m_n = 1.6749275e-27
M_sun = 1.98847e30
r0 = 1.25e-15

OUT_DIR = "outputs/tables"
os.makedirs(OUT_DIR, exist_ok=True)

lambdas = np.array([0.15, 0.18, 0.20, 0.22, 0.25, 0.28, 0.30, 0.32, 0.35])
mcrit_kg = np.sqrt((lambdas**3 * r0**3 * c**6) / (G**3 * m_n))
mcrit_msun = mcrit_kg / M_sun

df = pd.DataFrame({"lambda": lambdas, "Mcrit_kg": mcrit_kg, "Mcrit_Msun": mcrit_msun})
out = os.path.join(OUT_DIR, "scaling_threshold_table.csv")
df.to_csv(out, index=False)
print(df)
print(f"Saved: {out}")
