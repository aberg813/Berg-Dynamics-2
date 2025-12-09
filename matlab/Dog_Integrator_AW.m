function [SolutionToAlgebraicEquations, Output] = Dog_Integrator_AW
%DOG_INTEGRATOR_AW  Two-DOF launcher model derived from MotionGenesis equations.
%   This script keeps both generalized coordinates (phi for the handle and theta for
%   the sling) and enforces all assignment constraints: motor torque limit, total
%   length <= 0.5 m, theta within [-45 deg, 0 deg], release when the ball first
%   crosses the vertical axis, and realistic link masses/inertia from geometry.
%
%   Outputs:
%       SolutionToAlgebraicEquations.massMatrix(phi,theta)
%       SolutionToAlgebraicEquations.rhs(phi,theta,phidot,thetadot)
%           Anonymous functions that evaluate the symbolic EOM.
%       Output structure contains parameters, initial conditions, and release data.

%% Physical parameters (geometry ties directly to inertia and stiffness)
g        = 9.81;                 % m/s^2
E_cf     = 125e9;                % Pa (carbon fiber for link B)
rho_cf   = 1600;                 % kg/m^3
mBall    = 0.06;                 % kg
tau      = 0.7;                  % N*m (motor torque limit, applied to link A)

widthA   = 0.010;                % m
thickA   = 0.0015;               % m
rB       = 0.002;                % m (solid circular sling)
LA       = 0.15;                 % m
LB       = 0.35;                 % m  (LA + LB = 0.50 m)

areaA = widthA * thickA;
areaB = pi * rB^2;
mA    = rho_cf * areaA * LA;
mB    = rho_cf * areaB * LB;
IA    = mA * LA^2 / 12;
IB    = mB * LB^2 / 12;
JB    = 0.5 * pi * rB^4;
k     = E_cf * JB / LB;          % rotational spring stiffness

params = struct("g", g, "tau", tau, "LA", LA, "LB", LB, ...
    "mA", mA, "mB", mB, "mBall", mBall, "IA", IA, "IB", IB, "k", k);

%% Static equilibrium for theta (assignment item 5)
theta0 = findThetaEquilibrium(params);
if theta0 < deg2rad(-45) || theta0 > 0
    error("Equilibrium theta %.2f deg violates constraint. Adjust geometry.", rad2deg(theta0));
end

state0 = [0; theta0; 0; 0];      % [phi; theta; phiDot; thetaDot]
tSpan  = [0 1.5];

opts = odeset("RelTol", 1e-6, "AbsTol", 1e-8, "MaxStep", 1e-3, ...
    "Events", @(t,x) guardEvents(t, x, params));

[tSol, xSol, tEvt, xEvt, iEvt] = ode45(@(t,x) dynamics(t, x, params), tSpan, state0, opts);

%% Release detection and reporting
releaseIdx = find(iEvt == 1, 1, "last");
if isempty(releaseIdx)
    warning("No release detected before %.2f s. Check event settings.", tSpan(end));
    releaseState = [];
    releaseSpeed = NaN;
else
    releaseState = xEvt(releaseIdx, :);
    [~, vRel] = ballState(releaseState, params);
    releaseSpeed = norm(vRel);
    fprintf("Release at t = %.4f s, speed = %.3f m/s (phi = %.1f deg, theta = %.1f deg)\n", ...
        tEvt(releaseIdx), releaseSpeed, rad2deg(releaseState(1)), rad2deg(releaseState(2)));
end

plotResults(tSol, xSol, tEvt, iEvt, params, releaseIdx);

%% Package outputs for report reproducibility
SolutionToAlgebraicEquations = struct( ...
    "massMatrix", @(phi, theta) massMatrix(phi, theta, params), ...
    "rhs", @(phi, theta, phidot, thetadot) rhsVector(phi, theta, phidot, thetadot, params, params.tau));

Output = struct( ...
    "params", params, ...
    "theta0", theta0, ...
    "integrationTime", tSol, ...
    "stateHistory", xSol, ...
    "eventTimes", tEvt, ...
    "eventStates", xEvt, ...
    "eventIDs", iEvt, ...
    "releaseSpeed", releaseSpeed);

end

%--------------------------------------------------------------------------
function theta0 = findThetaEquilibrium(params)
% Solve k*theta + (0.5*mB + mBall)*g*LB*cos(theta) = 0 for theta in [-45,0] deg.
    g = params.g;
    LB = params.LB;
    mB = params.mB;
    mBall = params.mBall;
    k = params.k;
    balance = @(theta) -0.5*LB*g*mB*cos(theta) - LB*g*mBall*cos(theta) - k*theta;

    bracket = deg2rad([-60, -1]);
    if sign(balance(bracket(1))) == sign(balance(bracket(2)))
        error("Equilibrium bracket does not change sign. Adjust initial guess.");
    end
    theta0 = fzero(balance, bracket);
end

%--------------------------------------------------------------------------
function dx = dynamics(~, x, params)
% State derivative for [phi; theta; phiDot; thetaDot].
    phi      = x(1);
    theta    = x(2);
    phidot   = x(3);
    thetadot = x(4);

    M = massMatrix(phi, theta, params);
    rhs = rhsVector(phi, theta, phidot, thetadot, params, params.tau);
    acc = M \ rhs;

    dx = [phidot; thetadot; acc(1); acc(2)];
end

%--------------------------------------------------------------------------
function M = massMatrix(~, theta, params)
% Symmetric mass matrix from Lagrange derivation.
    LA = params.LA;
    LB = params.LB;
    mA = params.mA;
    mB = params.mB;
    mBall = params.mBall;
    IA = params.IA;
    IB = params.IB;

    cT = cos(theta);

    M11 = IA + IB + 0.25 * LA^2 * mA ...
        + 0.25 * mB * (4 * LA^2 + 4 * LA * LB * cT + LB^2) ...
        + mBall * (LA^2 + 2 * LA * LB * cT + LB^2);
    M12 = IB + 0.25 * LB * mB * (2 * LA * cT + LB) ...
        + LB * mBall * (LA * cT + LB);
    M22 = IB + 0.25 * LB^2 * mB + LB^2 * mBall;

    M = [M11, M12; M12, M22];
end

%--------------------------------------------------------------------------
function rhs = rhsVector(phi, theta, phidot, thetadot, params, tau)
% Right-hand side of M*qdd = rhs (gravity, Coriolis, spring, motor torque).
    LA = params.LA;
    LB = params.LB;
    mA = params.mA;
    mB = params.mB;
    mBall = params.mBall;
    g = params.g;
    k = params.k;

    cPhi = cos(phi);
    cPhiTheta = cos(phi + theta);
    sTheta = sin(theta);

    rhs1 = 0.5 * LA * LB * mB * thetadot * (2 * phidot + thetadot) * sTheta ...
        + LA * LB * mBall * thetadot * (2 * phidot + thetadot) * sTheta ...
        - 0.5 * LA * g * mA * cPhi ...
        - 0.5 * g * mB * (2 * LA * cPhi + LB * cPhiTheta) ...
        - g * mBall * (LA * cPhi + LB * cPhiTheta) ...
        + tau;

    rhs2 = -0.5 * LA * LB * mB * phidot^2 * sTheta ...
        - LA * LB * mBall * phidot^2 * sTheta ...
        - 0.5 * LB * g * mB * cPhiTheta ...
        - LB * g * mBall * cPhiTheta ...
        - k * theta;

    rhs = [rhs1; rhs2];
end

%--------------------------------------------------------------------------
function [value, isterminal, direction] = guardEvents(~, x, params)
% Event 1: ball crosses vertical axis (release, terminal)
% Event 2: theta < -45 deg (invalid)
% Event 3: theta > 0 deg (invalid)
% Event 4: phi < 0 deg (invalid, ensures CCW motion)
    phi   = x(1);
    theta = x(2);

    xBall = params.LA * cos(phi) + params.LB * cos(phi + theta);
    thetaLow = theta + deg2rad(45);
    thetaHigh = theta;
    phiBelow = phi + 1e-6;  % avoid triggering at t=0

    value = [xBall; thetaLow; thetaHigh; phiBelow];
    isterminal = [1; 1; 1; 1];
    direction = [-1; -1; 1; -1];
end

%--------------------------------------------------------------------------
function [pos, vel] = ballState(state, params)
% Position and velocity of the ball expressed in frame N.
    phi      = state(1);
    theta    = state(2);
    phidot   = state(3);
    thetadot = state(4);

    LA = params.LA;
    LB = params.LB;

    pos = [LA * cos(phi) + LB * cos(phi + theta);
           LA * sin(phi) + LB * sin(phi + theta)];

    vel = [-LA * phidot * sin(phi) - LB * (phidot + thetadot) * sin(phi + theta);
            LA * phidot * cos(phi) + LB * (phidot + thetadot) * cos(phi + theta)];
end

%--------------------------------------------------------------------------
function plotResults(t, x, tEvt, iEvt, params, releaseIdx)
% Angle/time histories plus ball speed profile.
    phi    = rad2deg(x(:, 1));
    theta  = rad2deg(x(:, 2));
    phidot = x(:, 3);
    thetadot = x(:, 4);

    releaseTime = [];
    if ~isempty(releaseIdx)
        releaseTime = tEvt(releaseIdx);
    end

    figure;
    subplot(2,1,1);
    plot(t, phi, "LineWidth", 1.5); hold on;
    plot(t, theta, "LineWidth", 1.5);
    if ~isempty(releaseTime)
        xline(releaseTime, "--k", "Release");
    end
    ylabel("Angle (deg)");
    legend("\phi", "\theta", "Location", "best");
    grid on;

    subplot(2,1,2);
    plot(t, phidot, "LineWidth", 1.5); hold on;
    plot(t, thetadot, "LineWidth", 1.5);
    if ~isempty(releaseTime)
        xline(releaseTime, "--k");
    end
    ylabel("Angular rate (rad/s)");
    xlabel("Time (s)");
    legend("\phi̇", "\thetȧ", "Location", "best");
    grid on;

    % Ball speed timeline
    speeds = zeros(size(t));
    for i = 1:numel(t)
        [~, v] = ballState(x(i, :), params);
        speeds(i) = norm(v);
    end
    figure;
    plot(t, speeds, "LineWidth", 1.5);
    if ~isempty(releaseTime)
        xline(releaseTime, "--k", "Release");
    end
    ylabel("Ball speed (m/s)");
    xlabel("Time (s)");
    grid on;
end
