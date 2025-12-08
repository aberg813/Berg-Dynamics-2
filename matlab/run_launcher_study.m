%% Dynamics of a two-link ball launcher
% This script reproduces the analysis documented in the project report.
clear; clc;

%% Step 1+2: Define geometry, materials, and assemble parameters
config.LA = 0.22;              % m
config.LB = 0.28;              % m (LA + LB <= 0.50 m)
config.widthA = 0.015;         % m
config.thicknessA = 0.002;     % m
config.rB = 8e-3;              % m
config.materialA = 'aluminum';
config.materialB = 'carbon_fiber';
config.torque = 0.7;           % N*m (motor limit)
config.g = 9.81;               % m/s^2
config.mBall = 0.06;           % kg

params = launcher_build_params(config);

paramTable = table({params.materialA.name; params.materialB.name}, ...
    [params.LA; params.LB], ...
    [params.mA; params.mB], ...
    [params.IA; params.IB], ...
    [params.crossA.moment; params.crossB.moment], ...
    [params.materialA.E; params.materialB.E], ...
    [NaN; params.k], ...
    'VariableNames', {'Material','Length_m','Mass_kg','Izz_kgm2','AreaMoment_m4','YoungsModulus_Pa','k_Nm_per_rad'}, ...
    'RowNames', {'Link A','Link B'});

disp('Parameter table (matches Step 2 deliverable):');
disp(paramTable);

%% Step 5: Static equilibrium for theta
balanceFun = @(th) launcher_torque_balance(th, params);
opt = optimoptions('fsolve','Display','none');
theta0 = fsolve(balanceFun, deg2rad(-5), opt);
theta0 = max(min(theta0, 0), deg2rad(-45));   % enforce physical window
fprintf('Static theta_0 = %.4f deg\n', rad2deg(theta0));

%% Step 6: Dynamic simulation with ode45
x0 = [0; 0; theta0; 0];
tSpan = [0 1.5];
odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-10,'Events', @(t,x) launcher_release_event(t,x,params));
[t, x, tRelease, xRelease] = ode45(@(t,x) launcher_dynamics(t,x,params), tSpan, x0, odeOpts);

if isempty(tRelease)
    error('Release event not detected. Increase torque or duration.');
end

[posRelease, velRelease] = launcher_ball_state(xRelease(end,:)', params);
releaseSpeed = norm(velRelease);
fprintf('Release achieved at t = %.3f s, v = %.3f m/s\n', tRelease(end), releaseSpeed);

% Kinematics history for Step 4
n = numel(t);
ballPos = zeros(n,2);
ballVel = zeros(n,2);
ballSpeed = zeros(n,1);
for k = 1:n
    [ballPos(k,:), ballVel(k,:)] = launcher_ball_state(x(k,:)', params);
    ballSpeed(k) = norm(ballVel(k,:));
end

%% Plotting (Step 6 deliverables)
figure;
subplot(2,1,1);
plot(t, rad2deg(x(:,1)), 'LineWidth', 1.5); hold on;
plot(t, rad2deg(x(:,3)), 'LineWidth', 1.5);
yline(rad2deg(theta0), '--k', 'Start');
xline(tRelease(end), '--k', 'Release');
ylabel('Angle (deg)'); legend('\phi','\theta'); grid on;
subplot(2,1,2);
plot(t, x(:,2), 'LineWidth', 1.5); hold on;
plot(t, x(:,4), 'LineWidth', 1.5);
xline(tRelease(end), '--k');
ylabel('Angular rate (rad/s)'); xlabel('Time (s)'); legend('\dot{\phi}','\dot{\theta}'); grid on;

grid on;

figure;
plot(ballPos(:,1), ballPos(:,2), 'LineWidth', 1.5); hold on;
scatter(ballPos(1,1), ballPos(1,2), 60, 'filled', 'g');
scatter(posRelease(1), posRelease(2), 60, 'filled', 'r');
xlabel('x (m)'); ylabel('y (m)'); axis equal; grid on;
title('Ball trajectory in frame N');

figure;
plot(t, ballSpeed, 'LineWidth', 1.5);
xline(tRelease(end), '--k');
ylabel('Speed (m/s)'); xlabel('Time (s)'); grid on;
title('Ball speed vs. time');

%% Step 7: Parametric sweep for release speed
LA_list = 0.20:0.02:0.30;
radius_list = [7e-3, 9e-3, 11e-3];
materials_B = {'aluminum','carbon_fiber','steel'};
maxRelease = -Inf;
records = {};

for LAi = LA_list
    LBi = 0.50 - LAi;
    if LBi <= 0.08, continue; end
    for rB = radius_list
        for mKey = materials_B
            cfg = config;
            cfg.LA = LAi; cfg.LB = LBi; cfg.rB = rB;
            cfg.materialB = mKey{1};
            params_i = launcher_build_params(cfg);
            balance_i = @(th) launcher_torque_balance(th, params_i);
            theta_guess = fsolve(balance_i, deg2rad(-10), opt);
            theta_guess = max(min(theta_guess, 0), deg2rad(-45));
            x0_i = [0; 0; theta_guess; 0];
            try
                [t_i, x_i, tRel_i, xRel_i] = ode45(@(t,x) launcher_dynamics(t,x,params_i), tSpan, x0_i, odeOpts);
                if isempty(tRel_i)
                    vRel = NaN; relTime = tSpan(end);
                else
                    [~, vel_i] = launcher_ball_state(xRel_i(end,:)', params_i);
                    vRel = norm(vel_i);
                    relTime = tRel_i(end);
                end
            catch
                vRel = NaN; relTime = NaN;
            end
            records(end+1,:) = {LAi, LBi, rB, params_i.materialB.name, rad2deg(theta_guess), vRel, relTime}; %#ok<AGROW>
            if ~isnan(vRel) && vRel > maxRelease
                maxRelease = vRel;
            end
        end
    end
end

sweepTable = cell2table(records, 'VariableNames', ...
    {'LA','LB','radius','materialB','theta0_deg','releaseSpeed','releaseTime'});

fprintf('Best release speed found in sweep: %.3f m/s\n', maxRelease);

%% Helper exports
writetable(paramTable, 'parameter_table.csv', 'WriteRowNames', true);
writetable(sweepTable, 'design_sweep_results.csv');
