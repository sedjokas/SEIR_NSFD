"""
SEIR NSFD vs. Standard Euler — Crisp and Fuzzy Simulations
Parameters adapted from Alqarni et al. (2023), amebiasis SEIR model.

Generates:
  AI2M4RI_PROCS_Template/seir_deterministic.pdf
  AI2M4RI_PROCS_Template/seir_fuzzy.pdf
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

OUT = "AI2M4RI_PROCS_Template"
os.makedirs(OUT, exist_ok=True)

# ── Crisp parameters (Alqarni et al. 2023) ────────────────────────────────────
LAM     = 0.02   # Λ  recruitment rate  (= μ → normalised DFE population N = 1)
MU      = 0.02   # μ  natural death rate
EPS     = 0.20   # ε  progression rate E → I
BETA    = 0.50   # β  transmission rate
D_CRIT  = 0.01   # d  disease-induced death rate
GAMMA   = 0.10   # γ  recovery rate

# Initial conditions (normalised, N = 1)
S0, E0, I0, R0 = 0.95, 0.03, 0.02, 0.00

T_END = 80   # simulation horizon

R0_val = BETA * LAM * EPS / (MU * (MU + EPS) * (MU + D_CRIT + GAMMA))
print(f"Basic reproduction number  R0 = {R0_val:.4f}")

LABELS = ["S  (Susceptible)", "E  (Exposed)", "I  (Infectious)", "R  (Recovered)"]


# ── Schemes ───────────────────────────────────────────────────────────────────
def nsfd_step(S, E, I, R, h, b=BETA, dv=D_CRIT, gv=GAMMA):
    Sn = (S + h * LAM)            / (1 + h * b * I   + h * MU)
    En = (E + h * b * Sn * I)     / (1 + h * (MU + EPS))
    In = (I + h * EPS * En)       / (1 + h * (MU + dv + gv))
    Rn = (R + h * gv * In)        / (1 + h * MU)
    return Sn, En, In, Rn


def euler_step(S, E, I, R, h, b=BETA, dv=D_CRIT, gv=GAMMA):
    dS = LAM - b * S * I - MU * S
    dE = b * S * I - (MU + EPS) * E
    dI = EPS * E - (MU + dv + gv) * I
    dR = gv * I - MU * R
    return S + h*dS, E + h*dE, I + h*dI, R + h*dR


def simulate(scheme, h, b=BETA, dv=D_CRIT, gv=GAMMA):
    steps = int(T_END / h)
    S, E, I, R = S0, E0, I0, R0
    traj = [(S, E, I, R)]
    for _ in range(steps):
        S, E, I, R = scheme(S, E, I, R, h, b, dv, gv)
        traj.append((S, E, I, R))
    t = np.linspace(0, T_END, steps + 1)
    return t, np.array(traj)


# ── Figure 1 : Crisp  NSFD vs. Euler ─────────────────────────────────────────
H_REF   = 0.01
H_LARGE = 2.0

t_ref,  tr  = simulate(nsfd_step,  H_REF)
t_nsfd, tn  = simulate(nsfd_step,  H_LARGE)
t_eul,  te  = simulate(euler_step, H_LARGE)

C_NSFD  = ["#1565C0", "#2E7D32", "#C62828", "#E65100"]
C_EULER = ["#90CAF9", "#A5D6A7", "#EF9A9A", "#FFCC80"]

fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True)
for k, ax in enumerate(axes.flat):
    ax.plot(t_ref,  tr[:, k], "k--",  lw=1.2, alpha=0.55, label=f"Reference (h={H_REF})")
    ax.plot(t_nsfd, tn[:, k], color=C_NSFD[k],  lw=2.2, label=f"NSFD   (h={H_LARGE})")
    ax.plot(t_eul,  te[:, k], color=C_EULER[k], lw=2.2, ls="-.", label=f"Euler  (h={H_LARGE})")
    ax.axhline(0, color="grey", lw=0.7, ls=":")
    ax.set_title(LABELS[k], fontsize=10, fontweight="bold")
    ax.set_ylabel("Population fraction", fontsize=8)
    ax.legend(fontsize=7.5)
    ax.set_ylim(-0.35, 1.12)

for ax in axes[1]:
    ax.set_xlabel("Time (days)", fontsize=9)

fig.suptitle(
    rf"SEIR amebiasis model — NSFD vs. Euler  ($h = {H_LARGE}$,  $\mathcal{{R}}_0 \approx {R0_val:.2f}$)",
    fontsize=11,
)
fig.tight_layout()
fig.savefig(f"{OUT}/seir_deterministic.pdf", bbox_inches="tight", dpi=150)
fig.savefig(f"{OUT}/seir_deterministic.png", bbox_inches="tight", dpi=150)
plt.close()
print("Figure 1 saved  →  seir_deterministic.pdf")


# ── Figure 2 : Fuzzy NSFD  (α-cut envelopes) ──────────────────────────────────
# Triangular fuzzy parameters  (a, m, b)
BETA_F  = (0.40, 0.50, 0.60)
D_F     = (0.008, 0.010, 0.012)
GAMMA_F = (0.080, 0.100, 0.120)

ALPHAS      = [0.0, 0.5, 1.0]
A_COLORS    = ["#B3E5FC", "#0288D1", "#01579B"]
H_FUZZY     = 0.5


def alpha_cut(tri, alpha):
    a, m, b = tri
    return a + alpha * (m - a), b - alpha * (b - m)


fig2, axes2 = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True)

# Reference crisp trajectory
t_c, tc = simulate(nsfd_step, H_FUZZY)

for idx, alpha in enumerate(ALPHAS):
    bl, bu = alpha_cut(BETA_F,  alpha)
    dl, du = alpha_cut(D_F,     alpha)
    gl, gu = alpha_cut(GAMMA_F, alpha)

    # Worst-case bounds for each compartment
    # Lower  I bound : low β, low d, high γ  → fewer infected
    # Upper  I bound : high β, high d, low γ → more infected
    _, tlo = simulate(nsfd_step, H_FUZZY, b=bl, dv=dl, gv=gu)
    _, thi = simulate(nsfd_step, H_FUZZY, b=bu, dv=du, gv=gl)
    t_a   = np.linspace(0, T_END, int(T_END / H_FUZZY) + 1)

    for k, ax in enumerate(axes2.flat):
        lo = np.minimum(tlo[:, k], thi[:, k])
        hi = np.maximum(tlo[:, k], thi[:, k])
        ax.fill_between(t_a, lo, hi,
                        alpha=0.30 + 0.18 * idx,
                        color=A_COLORS[idx])

for k, ax in enumerate(axes2.flat):
    ax.plot(t_c, tc[:, k], color="#01579B", lw=2, label=r"Crisp NSFD ($\alpha=1$)")
    ax.set_title(LABELS[k], fontsize=10, fontweight="bold")
    ax.set_ylabel("Population fraction", fontsize=8)
    ax.set_ylim(-0.02, 1.10)

for ax in axes2[1]:
    ax.set_xlabel("Time (days)", fontsize=9)

patches = [
    mpatches.Patch(color=A_COLORS[i], alpha=0.6 + 0.15 * i, label=rf"$\alpha = {a}$")
    for i, a in enumerate(ALPHAS)
]
crisp_line = plt.Line2D([0], [0], color="#01579B", lw=2, label=r"Crisp NSFD ($\alpha=1$)")
axes2.flat[0].legend(handles=patches + [crisp_line], fontsize=8, loc="upper right")

fig2.suptitle(
    r"Fuzzy NSFD — $\alpha$-cut envelopes  ($h=0.5$)"
    "\n"
    r"$\tilde{\beta}=(0.40,\,0.50,\,0.60)$,  "
    r"$\tilde{d}=(0.008,\,0.010,\,0.012)$,  "
    r"$\tilde{\gamma}=(0.08,\,0.10,\,0.12)$",
    fontsize=10,
)
fig2.tight_layout()
fig2.savefig(f"{OUT}/seir_fuzzy.pdf", bbox_inches="tight", dpi=150)
fig2.savefig(f"{OUT}/seir_fuzzy.png", bbox_inches="tight", dpi=150)
plt.close()
print("Figure 2 saved  →  seir_fuzzy.pdf")
print("Done.")
