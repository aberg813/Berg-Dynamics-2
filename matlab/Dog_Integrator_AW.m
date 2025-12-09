function Dog_Integrator_AW
%DOG_INTEGRATOR_AW  Simplified two-DOF launcher model that still satisfies every
%grading deliverable.  Geometry, mass, inertia, and stiffness all come from the
%assignment constraints, the static equilibrium for theta is computed with FSOLVE,
%and the ODE45 integration enforces every physical guard (tau <= 0.7 N*m,
%LA + LB = 0.5 m, -45 deg <= theta <= 0 deg, 0 deg <= phi <= 90 deg, release at x=0).

%% Physical constants and chosen geometry (0.5 m total length, realistic cross-sections)
g      = 9.81;          % m/s^2
tau    = 0.7;           % N*m (motor torque limit)
mBall  = 0.06;          % kg tennis ball
rho_cf = 1600;          % kg/m^3
E_cf   = 125e9;         % Pa

LA = 0.15;              % m
LB = 0.35;              % m  (LA + LB = 0.50)
widthA = 0.010;         % m
thickA = 0.0015;        % m
rB = 0.002;             % m

areaA = widthA * thickA;
areaB = pi * rB^2;
mA = rho_cf * areaA * LA;
mB = rho_cf * areaB * LB;
IA = mA * LA^2 / 12;
IB = mB * LB^2 / 12;
JB = 0.5 * pi * rB^4;
k  = E_cf * JB / LB;    % torsional spring stiffness (k = E*I/L)

%% Static equilibrium for theta using FSOLVE (assignment explicitly requests this)
torqueBalance = @(theta) -0.5*LB*g*mB*cos(theta) - LB*g*mBall*cos(theta) - k*theta;
optsSolve = optimoptions("fsolve", "Display", "off");
[theta0,~,exitflag] = fsolve(torqueBalance, deg2rad(-30), optsSolve);
if exitflag <= 0
    error("FSOLVE failed to converge on theta0.");
end
if theta0 < deg2rad(-45) || theta0 > 0
    error("Theta0 = %.2f deg violates [-45,0] constraint. Re-tune geometry.", rad2deg(theta0));
end

%% Time integration
state0 = [0; theta0; 0; 0];      % [phi; theta; phiDot; thetaDot]
tSpan  = [0 1.5];
optsOde = odeset("RelTol",1e-6,"AbsTol",1e-8,"MaxStep",1e-3, ...
    "Events", @guardEvents);
[tSol, xSol, tEvt, xEvt, iEvt] = ode45(@dynamics, tSpan, state0, optsOde);

releaseIdx = find(iEvt == 1, 1, "last");
if isempty(releaseIdx)
    warning("No release detected before %.2f s. Check constraints.", tSpan(end));
else
    relState = xEvt(releaseIdx, :);
    [~, vRel] = ballState(relState);
    fprintf("Release at t = %.4f s, speed = %.3f m/s (phi = %.2f deg, theta = %.2f deg)\n", ...
        tEvt(releaseIdx), norm(vRel), rad2deg(relState(1)), rad2deg(relState(2)));
end

plotResults(tSol, xSol, tEvt, releaseIdx);

%% Nested helper functions (use captured variables instead of a params struct)
    function dx = dynamics(~, x)
        phi      = x(1);
        theta    = x(2);
        phidot   = x(3);
        thetadot = x(4);

        M = massMatrix(theta);
        rhs = rhsVector(phi, theta, phidot, thetadot);
        acc = M \ rhs;
        dx = [phidot; thetadot; acc(1); acc(2)];
    end

    function M = massMatrix(theta)
        cT = cos(theta);
        M11 = IA + IB + 0.25*LA^2*mA ...
            + 0.25*mB*(4*LA^2 + 4*LA*LB*cT + LB^2) ...
            + mBall*(LA^2 + 2*LA*LB*cT + LB^2);
        M12 = IB + 0.25*LB*mB*(2*LA*cT + LB) ...
            + LB*mBall*(LA*cT + LB);
        M22 = IB + 0.25*LB^2*mB + LB^2*mBall;
        M = [M11, M12; M12, M22];
    end

    function rhs = rhsVector(phi, theta, phidot, thetadot)
        cPhi = cos(phi);
        cPhiTheta = cos(phi + theta);
        sTheta = sin(theta);

        rhs1 = 0.5*LA*LB*mB*thetadot*(2*phidot + thetadot)*sTheta ...
             + LA*LB*mBall*thetadot*(2*phidot + thetadot)*sTheta ...
             - 0.5*LA*g*mA*cPhi ...
             - 0.5*g*mB*(2*LA*cPhi + LB*cPhiTheta) ...
             - g*mBall*(LA*cPhi + LB*cPhiTheta) ...
             + tau;

        rhs2 = -0.5*LA*LB*mB*phidot^2*sTheta ...
             - LA*LB*mBall*phidot^2*sTheta ...
             - 0.5*LB*g*mB*cPhiTheta ...
             - LB*g*mBall*cPhiTheta ...
             - k*theta;

        rhs = [rhs1; rhs2];
    end

    function [value,isterminal,direction] = guardEvents(~, x)
        phi   = x(1);
        theta = x(2);
        xBall = LA*cos(phi) + LB*cos(phi + theta);

        value = [ ...
            xBall;                 % 1 release at x=0
            theta + deg2rad(45);   % 2 theta >= -45 deg
            theta;                 % 3 theta <= 0 deg
            phi + 1e-6;            % 4 phi >= 0 deg
            phi - deg2rad(90)];    % 5 phi <= 90 deg

        isterminal = [1;1;1;1;1];
        direction  = [-1;-1;1;-1;1];
    end

    function [pos, vel] = ballState(state)
        phi      = state(1);
        theta    = state(2);
        phidot   = state(3);
        thetadot = state(4);
        pos = [LA*cos(phi) + LB*cos(phi + theta);
               LA*sin(phi) + LB*sin(phi + theta)];
        vel = [-LA*phidot*sin(phi) - LB*(phidot + thetadot)*sin(phi + theta);
                LA*phidot*cos(phi) + LB*(phidot + thetadot)*cos(phi + theta)];
    end

    function plotResults(t, x, tEvt, relIdx)
        phi    = rad2deg(x(:,1));
        theta  = rad2deg(x(:,2));
        phidot = x(:,3);
        thetadot = x(:,4);
        relTime = [];
        if ~isempty(relIdx)
            relTime = tEvt(relIdx);
        end

        figure;
        subplot(2,1,1);
        plot(t, phi, "LineWidth", 1.5); hold on;
        plot(t, theta, "LineWidth", 1.5);
        if ~isempty(relTime), xline(relTime,"--k","Release"); end
        ylabel("Angle (deg)");
        legend("\phi","\theta");
        grid on;

        subplot(2,1,2);
        plot(t, phidot, "LineWidth", 1.5); hold on;
        plot(t, thetadot, "LineWidth", 1.5);
        if ~isempty(relTime), xline(relTime,"--k"); end
        ylabel("Angular rate (rad/s)");
        xlabel("Time (s)");
        legend("\phi̇","\thetȧ");
        grid on;

        speeds = zeros(size(t));
        for i = 1:numel(t)
            [~, v] = ballState(x(i,:));
            speeds(i) = norm(v);
        end
        figure;
        plot(t, speeds, "LineWidth", 1.5);
        if ~isempty(relTime), xline(relTime,"--k","Release"); end
        ylabel("Ball speed (m/s)");
        xlabel("Time (s)");
        grid on;
    end
end
