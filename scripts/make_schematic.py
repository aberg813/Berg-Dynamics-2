import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc, Rectangle, Circle

DATA_PATH = Path("data/base_summary.json")
FIG_PATH = Path("figures/schematic.png")

with open(DATA_PATH, "r", encoding="utf-8") as f:
    summary = json.load(f)

LA = summary["geometry"]["LA"]
LB = summary["geometry"]["LB"]
radius_B = 0.008  # match GEOMETRY configuration
width_A = 0.015
thickness_A = 0.002

phi_demo = np.deg2rad(35)
theta_demo = np.deg2rad(-12)

base = np.array([0.0, 0.0])
joint = np.array([LA * np.cos(phi_demo), LA * np.sin(phi_demo)])
com_A = base + 0.5 * (joint - base)
linkB_dir = np.array([np.cos(phi_demo + theta_demo), np.sin(phi_demo + theta_demo)])
com_B = joint + 0.5 * LB * linkB_dir
ball = joint + LB * linkB_dir

fig = plt.figure(figsize=(8, 5))
gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.3)
ax = fig.add_subplot(gs[0, 0])
ax_cross = fig.add_subplot(gs[0, 1])

ax.plot([base[0], joint[0]], [base[1], joint[1]], color="#1f77b4", linewidth=6, label="Link A")
ax.plot([joint[0], ball[0]], [joint[1], ball[1]], color="#ff7f0e", linewidth=6, label="Link B")
ax.scatter([base[0], joint[0], com_A[0], com_B[0], ball[0]],
           [base[1], joint[1], com_A[1], com_B[1], ball[1]],
           color=["black", "black", "#1f77b4", "#ff7f0e", "gold"], zorder=5)
ax.text(com_A[0], com_A[1] + 0.02, "C_A", color="#1f77b4")
ax.text(com_B[0], com_B[1] - 0.025, "C_B", color="#ff7f0e")
ax.text(ball[0] + 0.01, ball[1], "Ball", color="goldenrod")

# axes
i_vec = np.array([0.1, 0])
j_vec = np.array([0, 0.1])
ax.arrow(-0.02, -0.02, i_vec[0], i_vec[1], head_width=0.01, color="k")
ax.arrow(-0.02, -0.02, j_vec[0], j_vec[1], head_width=0.01, color="k")
ax.text(0.09, -0.015, r"$\hat{i}_N$")
ax.text(-0.015, 0.09, r"$\hat{j}_N$")

# angle phi arc
arc_phi = Arc(base, 0.25, 0.25, theta1=0, theta2=np.degrees(phi_demo), color="k")
ax.add_patch(arc_phi)
ax.text(0.11, 0.03, r"$\phi$")
# theta arc about joint
arc_theta = Arc(joint, 0.18, 0.18, angle=np.degrees(phi_demo),
                theta1=0, theta2=np.degrees(theta_demo), color="k")
ax.add_patch(arc_theta)
ax.text(joint[0] + 0.05, joint[1] - 0.01, r"$\theta$")

# annotate pivot and spring
ax.text(base[0] - 0.015, base[1] - 0.02, "motor/torque τ")
ax.annotate("rotational spring k", xy=joint, xytext=(joint[0] + 0.05, joint[1] - 0.1),
            arrowprops=dict(arrowstyle="->"))

ax.set_aspect('equal')
ax.set_xlim(-0.1, 0.4)
ax.set_ylim(-0.05, 0.35)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.legend(loc="upper left")
ax.set_title("Two-link launcher schematic")
ax.grid(True, linestyle=":", alpha=0.4)

# cross-section subplot
ax_cross.add_patch(Circle((0.0, 0.0), radius_B, facecolor="#ff7f0e", alpha=0.6))
ax_cross.set_aspect('equal')
ax_cross.set_xlim(-0.012, 0.012)
ax_cross.set_ylim(-0.012, 0.012)
ax_cross.set_xticks([])
ax_cross.set_yticks([])
ax_cross.set_title("Link B cross-section")
ax_cross.annotate(f"r = {radius_B*1000:.0f} mm", xy=(0.0, radius_B), xytext=(0.004, 0.010),
                  arrowprops=dict(arrowstyle="->"))
ax_cross.text(-0.0115, -0.0115, "solid CF rod", fontsize=9)

fig.suptitle("Ball launcher geometry and symbols", fontsize=14)
fig.tight_layout()
fig.savefig(FIG_PATH, dpi=300)
plt.close(fig)
