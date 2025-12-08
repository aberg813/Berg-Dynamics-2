import json
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt


@dataclass
class Material:
    name: str
    density: float  # kg/m^3
    youngs_modulus: float  # Pa


@dataclass
class Geometry:
    kind: str
    params: Dict[str, float]


MATERIALS = {
    "aluminum": Material("Aluminum 6061", density=2700.0, youngs_modulus=69e9),
    "steel": Material("AISI 1080 Steel", density=7850.0, youngs_modulus=200e9),
    "carbon_fiber": Material("Unidirectional Carbon Fiber", density=1600.0, youngs_modulus=125e9),
}

GEOMETRY = {
    "A": Geometry("rect", {"width": 0.015, "thickness": 0.002}),
    "B": Geometry("circle", {"radius": 0.008}),
}

g = 9.81
TAU_MAX = 0.7  # N*m constant motor torque applied to link A

SYMBOLIC_CACHE = {}


def build_symbolic_model():
    """Construct and cache the symbolic equations of motion and kinematic maps."""
    if SYMBOLIC_CACHE:
        return SYMBOLIC_CACHE

    t = sp.symbols("t")
    phi_fun = sp.Function("phi")(t)
    theta_fun = sp.Function("theta")(t)
    phidot_fun = sp.diff(phi_fun, t)
    thetadot_fun = sp.diff(theta_fun, t)
    phiddot_fun = sp.diff(phi_fun, t, 2)
    thetaddot_fun = sp.diff(theta_fun, t, 2)

    mA, mB, mBall, LA, LB = sp.symbols("mA mB mBall LA LB", positive=True)
    IA, IB = sp.symbols("IA IB", positive=True)
    k = sp.symbols("k", positive=True)
    tau = sp.symbols("tau")
    grav = sp.symbols("g", positive=True)

    cA = LA / 2
    cB = LB / 2

    rA = sp.Matrix([
        cA * sp.cos(phi_fun),
        cA * sp.sin(phi_fun)
    ])
    rJoint = sp.Matrix([
        LA * sp.cos(phi_fun),
        LA * sp.sin(phi_fun)
    ])
    rB = rJoint + sp.Matrix([
        cB * sp.cos(phi_fun + theta_fun),
        cB * sp.sin(phi_fun + theta_fun)
    ])
    rBall = rJoint + sp.Matrix([
        LB * sp.cos(phi_fun + theta_fun),
        LB * sp.sin(phi_fun + theta_fun)
    ])

    vA = sp.diff(rA, t)
    vB = sp.diff(rB, t)
    vBall = sp.diff(rBall, t)

    omegaA = phidot_fun
    omegaB = phidot_fun + thetadot_fun

    T = (
        0.5 * mA * (vA.dot(vA)) + 0.5 * IA * omegaA ** 2 +
        0.5 * mB * (vB.dot(vB)) + 0.5 * IB * omegaB ** 2 +
        0.5 * mBall * (vBall.dot(vBall))
    )
    V = mA * grav * rA[1] + mB * grav * rB[1] + mBall * grav * rBall[1] + 0.5 * k * theta_fun ** 2
    L = T - V

    qs = [phi_fun, theta_fun]
    Qs = [tau, 0]
    eqs = []
    for qi, Qi in zip(qs, Qs):
        eq = sp.diff(sp.diff(L, sp.diff(qi, t)), t) - sp.diff(L, qi) - Qi
        eqs.append(sp.simplify(eq))

    phi, theta, phidot, thetadot, phidd, thetadd = sp.symbols("phi theta phidot thetadot phidd thetadd")
    subs = {
        phi_fun: phi,
        theta_fun: theta,
        phidot_fun: phidot,
        thetadot_fun: thetadot,
        phiddot_fun: phidd,
        thetaddot_fun: thetadd,
    }
    eqs_plain = [eq.subs(subs) for eq in eqs]
    M, forcing = sp.linear_eq_to_matrix(eqs_plain, [phidd, thetadd])

    M_func = sp.lambdify((phi, theta, mA, mB, mBall, LA, LB, IA, IB), M, modules="numpy")
    forcing_func = sp.lambdify((phi, theta, phidot, thetadot, mA, mB, mBall, LA, LB, grav, k, tau), forcing, modules="numpy")

    x_ball = rBall[0].subs({phi_fun: phi, theta_fun: theta})
    y_ball = rBall[1].subs({phi_fun: phi, theta_fun: theta})
    vx_ball = vBall[0].subs({phi_fun: phi, theta_fun: theta, phidot_fun: phidot, thetadot_fun: thetadot})
    vy_ball = vBall[1].subs({phi_fun: phi, theta_fun: theta, phidot_fun: phidot, thetadot_fun: thetadot})

    x_ball_func = sp.lambdify((phi, theta, LA, LB), x_ball, modules="numpy")
    y_ball_func = sp.lambdify((phi, theta, LA, LB), y_ball, modules="numpy")
    vx_ball_func = sp.lambdify((phi, theta, phidot, thetadot, LA, LB), vx_ball, modules="numpy")
    vy_ball_func = sp.lambdify((phi, theta, phidot, thetadot, LA, LB), vy_ball, modules="numpy")

    dV_dtheta = sp.diff(V, theta_fun).subs({phi_fun: phi, theta_fun: theta})
    torque_balance_func = sp.lambdify((theta, phi, mA, mB, mBall, LA, LB, grav, k), dV_dtheta, modules="numpy")

    SYMBOLIC_CACHE.update(
        M_func=M_func,
        forcing_func=forcing_func,
        x_ball=x_ball_func,
        y_ball=y_ball_func,
        vx_ball=vx_ball_func,
        vy_ball=vy_ball_func,
        torque_balance=torque_balance_func,
    )
    return SYMBOLIC_CACHE


def compute_area_moment(geometry: Geometry) -> float:
    if geometry.kind == "rect":
        b = geometry.params["width"]
        h = geometry.params["thickness"]
        return (b * h / 12.0) * (b ** 2 + h ** 2)
    if geometry.kind == "circle":
        r = geometry.params["radius"]
        return 0.5 * np.pi * r ** 4
    raise ValueError("Unsupported geometry")


def compute_cross_section_area(geometry: Geometry) -> float:
    if geometry.kind == "rect":
        return geometry.params["width"] * geometry.params["thickness"]
    if geometry.kind == "circle":
        r = geometry.params["radius"]
        return np.pi * r ** 2
    raise ValueError("Unsupported geometry")


def compute_link_properties(material_key: str, length: float, geometry: Geometry) -> Dict[str, float]:
    material = MATERIALS[material_key]
    area = compute_cross_section_area(geometry)
    volume = area * length
    mass = material.density * volume
    Izz = (mass * length ** 2) / 12.0
    area_moment = compute_area_moment(geometry)
    k = material.youngs_modulus * area_moment / length
    return {
        "material": material.name,
        "density": material.density,
        "E": material.youngs_modulus,
        "area": area,
        "volume": volume,
        "mass": mass,
        "Izz": Izz,
        "area_moment": area_moment,
        "k": k,
    }


def find_theta_equilibrium(params: Dict[str, float]) -> float:
    funcs = build_symbolic_model()
    balance = funcs["torque_balance"]

    def f(theta):
        return balance(theta, 0.0, params["mA"], params["mB"], params["mBall"], params["LA"], params["LB"], g, params["k"])

    result = root_scalar(f, bracket=(-0.9, 0.1), method="brentq")
    if not result.converged:
        raise RuntimeError("Failed to find equilibrium theta")
    return result.root


def simulate_release(params: Dict[str, float], theta0: float, t_final: float = 1.5):
    funcs = build_symbolic_model()
    M_func = funcs["M_func"]
    forcing_func = funcs["forcing_func"]
    x_ball = funcs["x_ball"]
    y_ball = funcs["y_ball"]
    vx_ball = funcs["vx_ball"]
    vy_ball = funcs["vy_ball"]

    def dynamics(t, state):
        phi, phidot, theta, thetadot = state
        M = np.array(M_func(phi, theta, params["mA"], params["mB"], params["mBall"], params["LA"], params["LB"], params["IA"], params["IB"]))
        forcing = np.array(forcing_func(phi, theta, phidot, thetadot, params["mA"], params["mB"], params["mBall"], params["LA"], params["LB"], g, params["k"], TAU_MAX)).reshape(2)
        acc = np.linalg.solve(M, forcing)
        return [phidot, acc[0], thetadot, acc[1]]

    def release_event(t, state):
        phi, _, theta, _ = state
        return x_ball(phi, theta, params["LA"], params["LB"])

    release_event.direction = -1
    release_event.terminal = True

    state0 = [0.0, 0.0, theta0, 0.0]
    sol = solve_ivp(dynamics, [0, t_final], state0, max_step=1e-3, events=release_event, rtol=1e-6, atol=1e-9)

    if len(sol.t_events[0]) == 0:
        raise RuntimeError("Release event not detected within simulation time")

    t_release = sol.t_events[0][0]
    state_release = sol.y_events[0][0]
    phi_r, phidot_r, theta_r, thetadot_r = state_release

    vx = vx_ball(phi_r, theta_r, phidot_r, thetadot_r, params["LA"], params["LB"])
    vy = vy_ball(phi_r, theta_r, phidot_r, thetadot_r, params["LA"], params["LB"])
    speed = float(np.sqrt(vx ** 2 + vy ** 2))

    traj = {
        "time": sol.t.tolist(),
        "phi": sol.y[0].tolist(),
        "phidot": sol.y[1].tolist(),
        "theta": sol.y[2].tolist(),
        "thetadot": sol.y[3].tolist(),
        "x_ball": [float(x_ball(phi, theta, params["LA"], params["LB"])) for phi, theta in zip(sol.y[0], sol.y[2])],
        "y_ball": [float(y_ball(phi, theta, params["LA"], params["LB"])) for phi, theta in zip(sol.y[0], sol.y[2])],
    }

    release_state = {
        "time": t_release,
        "phi": phi_r,
        "theta": theta_r,
        "phidot": phidot_r,
        "thetadot": thetadot_r,
        "speed": speed,
        "vx": float(vx),
        "vy": float(vy),
    }
    return release_state, traj


def assemble_parameters(LA: float, LB: float, geom_B: Geometry, material_A: str = "aluminum", material_B: str = "carbon_fiber"):
    linkA = compute_link_properties(material_A, LA, GEOMETRY["A"])
    linkB = compute_link_properties(material_B, LB, geom_B)
    params = {
        "LA": LA,
        "LB": LB,
        "mA": linkA["mass"],
        "mB": linkB["mass"],
        "IA": linkA["Izz"],
        "IB": linkB["Izz"],
        "k": linkB["k"],
        "mBall": 0.06,
        "linkA": linkA,
        "linkB": linkB,
    }
    return params


def export_time_series(traj: Dict[str, List[float]], base_path: Path):
    df = pd.DataFrame(traj)
    df.to_csv(base_path / "time_histories.csv", index=False)


def plot_timeseries(traj: Dict[str, List[float]], release_state: Dict[str, float], params: Dict[str, float], fig_dir: Path):
    t = traj["time"]
    phi = np.array(traj["phi"]) * 180 / np.pi
    theta = np.array(traj["theta"]) * 180 / np.pi
    phidot = np.array(traj["phidot"])
    thetadot = np.array(traj["thetadot"])
    funcs = build_symbolic_model()
    vx_ball_func = funcs["vx_ball"]
    vy_ball_func = funcs["vy_ball"]

    fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    ax[0].plot(t, phi, label="ϕ (deg)")
    ax[0].plot(t, theta, label="θ (deg)")
    ax[0].axvline(release_state["time"], color="k", linestyle="--", label="Release")
    ax[0].set_ylabel("Angle (deg)")
    ax[0].legend()
    ax[0].grid(True)

    ax[1].plot(t, phidot, label="ϕ̇ (rad/s)")
    ax[1].plot(t, thetadot, label="θ̇ (rad/s)")
    ax[1].axvline(release_state["time"], color="k", linestyle="--")
    ax[1].set_ylabel("Angular rate (rad/s)")
    ax[1].set_xlabel("Time (s)")
    ax[1].grid(True)
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "angles_rates.png", dpi=300)
    plt.close(fig)

    x = traj["x_ball"]
    y = traj["y_ball"]
    y_release = funcs["y_ball"](release_state["phi"], release_state["theta"], params["LA"], params["LB"])
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(x, y, label="Ball path")
    ax.scatter(x[0], y[0], color="green", label="Start")
    ax.scatter(0.0, float(y_release), color="red", label="Release")
    ax.axvline(0, color="k", linestyle=":")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Ball trajectory in frame N")
    ax.grid(True)
    ax.axis("equal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "trajectory.png", dpi=300)
    plt.close(fig)

    speeds = []
    for phi_i, theta_i, phidot_i, thetadot_i in zip(traj["phi"], traj["theta"], traj["phidot"], traj["thetadot"]):
        vx = vx_ball_func(phi_i, theta_i, phidot_i, thetadot_i, params["LA"], params["LB"])
        vy = vy_ball_func(phi_i, theta_i, phidot_i, thetadot_i, params["LA"], params["LB"])
        speeds.append(np.sqrt(vx ** 2 + vy ** 2))
    speeds = np.array(speeds)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(t, speeds)
    ax.axvline(release_state["time"], color="k", linestyle="--")
    ax.set_ylabel("Ball speed (m/s)")
    ax.set_xlabel("Time (s)")
    ax.set_title("Ball speed profile")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(fig_dir / "generalized_speed.png", dpi=300)
    plt.close(fig)


def plot_torque_balance(theta_vals: np.ndarray, balance_vals: np.ndarray, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(theta_vals * 180 / np.pi, balance_vals)
    ax.axhline(0, color="k", linestyle=":")
    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("Torque balance f(θ)")
    ax.set_title("Static equilibrium search")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(fig_dir / "torque_balance.png", dpi=300)
    plt.close(fig)


def sweep_design_space(LA_values, radius_values, materials_B) -> pd.DataFrame:
    records = []
    funcs = build_symbolic_model()
    x_ball = funcs["x_ball"]

    for LA in LA_values:
        LB = 0.5 - LA
        if LB <= 0.05:
            continue
        for radius in radius_values:
            geom_B = Geometry("circle", {"radius": radius})
            for material_key in materials_B:
                params = assemble_parameters(LA, LB, geom_B, material_B=material_key)
                try:
                    theta0 = find_theta_equilibrium(params)
                    release_state, _ = simulate_release(params, theta0)
                    record = {
                        "LA": LA,
                        "LB": LB,
                        "radius": radius,
                        "material_B": MATERIALS[material_key].name,
                        "theta0_deg": np.degrees(theta0),
                        "release_speed": release_state["speed"],
                        "release_time": release_state["time"],
                    }
                    records.append(record)
                except Exception as exc:  # noqa: PERF203
                    record = {
                        "LA": LA,
                        "LB": LB,
                        "radius": radius,
                        "material_B": MATERIALS[material_key].name,
                        "theta0_deg": np.nan,
                        "release_speed": np.nan,
                        "release_time": np.nan,
                        "note": str(exc),
                    }
                    records.append(record)
    df = pd.DataFrame(records)
    return df


def plot_optimization_results(df: pd.DataFrame, fig_dir: Path):
    df_valid = df.dropna(subset=["release_speed"])
    fig, ax = plt.subplots(figsize=(7, 4))
    markers = ["o", "s", "^", "D", "P"]
    scatter_ref = None
    for idx, material in enumerate(df_valid["material_B"].unique()):
        marker = markers[idx % len(markers)]
        mask = df_valid["material_B"] == material
        scatter = ax.scatter(df_valid.loc[mask, "LA"], df_valid.loc[mask, "release_speed"],
                             c=df_valid.loc[mask, "radius"], cmap="viridis", s=80,
                             marker=marker, edgecolors="k", label=material, alpha=0.85)
        if scatter_ref is None:
            scatter_ref = scatter
    ax.set_xlabel("L_A (m)")
    ax.set_ylabel("Release speed (m/s)")
    ax.set_title("Design sweep results")
    if scatter_ref is not None:
        fig.colorbar(scatter_ref, label="Link B radius (m)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(fig_dir / "design_sweep.png", dpi=300)
    plt.close(fig)


def main():
    base_geom_B = GEOMETRY["B"]
    base_params = assemble_parameters(0.22, 0.28, base_geom_B)

    theta0 = find_theta_equilibrium(base_params)
    release_state, traj = simulate_release(base_params, theta0)

    out_dir = Path("data")
    fig_dir = Path("figures")
    out_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(exist_ok=True)

    export_time_series(traj, out_dir)
    plot_timeseries(traj, release_state, base_params, fig_dir)

    funcs = build_symbolic_model()
    balance = funcs["torque_balance"]
    theta_vals = np.linspace(-np.pi / 2, 0.2, 200)
    balance_vals = [balance(th, 0.0, base_params["mA"], base_params["mB"], base_params["mBall"], base_params["LA"], base_params["LB"], g, base_params["k"]) for th in theta_vals]
    plot_torque_balance(theta_vals, np.array(balance_vals, dtype=float), fig_dir)

    summary = {
        "theta0_deg": float(np.degrees(theta0)),
        "release_speed": release_state["speed"],
        "release_time": release_state["time"],
        "geometry": {
            "LA": base_params["LA"],
            "LB": base_params["LB"],
        },
        "masses": {
            "mA": base_params["mA"],
            "mB": base_params["mB"],
            "mBall": base_params["mBall"],
        },
        "release_angles_deg": {
            "phi": float(np.degrees(release_state["phi"])),
            "theta": float(np.degrees(release_state["theta"])),
        },
        "release_rates": {
            "phidot": release_state["phidot"],
            "thetadot": release_state["thetadot"],
        },
        "linkA": base_params["linkA"],
        "linkB": base_params["linkB"],
    }

    with open(out_dir / "base_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    df_params = pd.DataFrame([
        {
            "Link": "A",
            "Material": base_params["linkA"]["material"],
            "Length_m": base_params["LA"],
            "Mass_kg": base_params["mA"],
            "Izz_kgm2": base_params["IA"],
            "AreaMoment_m4": base_params["linkA"]["area_moment"],
            "E_Pa": base_params["linkA"]["E"],
            "k_Nm_rad": np.nan,
        },
        {
            "Link": "B",
            "Material": base_params["linkB"]["material"],
            "Length_m": base_params["LB"],
            "Mass_kg": base_params["mB"],
            "Izz_kgm2": base_params["IB"],
            "AreaMoment_m4": base_params["linkB"]["area_moment"],
            "E_Pa": base_params["linkB"]["E"],
            "k_Nm_rad": base_params["k"],
        }
    ])
    df_params.to_csv(out_dir / "parameter_table.csv", index=False)

    LA_values = np.linspace(0.22, 0.30, 5)
    radius_values = np.array([0.009, 0.011, 0.013])
    materials_B = ["aluminum", "carbon_fiber", "steel"]
    df_sweep = sweep_design_space(LA_values, radius_values, materials_B)
    df_sweep.to_csv(out_dir / "design_sweep.csv", index=False)
    plot_optimization_results(df_sweep, fig_dir)

    best_row = df_sweep.sort_values("release_speed", ascending=False).dropna(subset=["release_speed"]).iloc[0]
    with open(out_dir / "best_design.json", "w", encoding="utf-8") as f:
        json.dump(best_row.to_dict(), f, indent=2)


if __name__ == "__main__":
    main()
