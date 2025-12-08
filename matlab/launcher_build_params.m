function params = launcher_build_params(config)
%LAUNCHER_BUILD_PARAMS  Assemble physical parameters for the launcher.
%
% Inputs:
%   config.LA, config.LB  - link lengths (m)
%   config.widthA, config.thicknessA - link A rectangular cross-section (m)
%   config.rB             - link B circular radius (m)
%   config.materialA/B    - material strings: 'aluminum', 'steel', 'carbon_fiber'
%   config.torque         - applied motor torque (N*m)
%   config.mBall          - ball mass (kg)
%   config.g              - gravitational acceleration (m/s^2)
%
% Outputs:
%   params struct with masses, inertias, stiffness, and metadata.

mats = launcher_material_library();
matA = mats.(lower(config.materialA));
matB = mats.(lower(config.materialB));

areaA = config.widthA * config.thicknessA;
volA = areaA * config.LA;
mA = matA.density * volA;
IA = (mA * config.LA^2) / 12;
JA = (areaA / 12) * (config.widthA^2 + config.thicknessA^2); % polar approx

areaB = pi * config.rB^2;
volB = areaB * config.LB;
mB = matB.density * volB;
IB = (mB * config.LB^2) / 12;
JB = 0.5 * pi * config.rB^4;
k = matB.E * JB / config.LB;

params = struct();
params.LA = config.LA;
params.LB = config.LB;
params.mA = mA;
params.mB = mB;
params.mBall = config.mBall;
params.IA = IA;
params.IB = IB;
params.crossA = struct('area', areaA, 'moment', JA);
params.crossB = struct('area', areaB, 'moment', JB, 'radius', config.rB);
params.materialA = matA;
params.materialB = matB;
params.k = k;
params.tau = config.torque;
params.g = config.g;
params.geometry = config;
end

function mats = launcher_material_library()
% Material properties (density kg/m^3, Young's modulus Pa)
mat.aluminum = struct('name', 'Aluminum 6061', 'density', 2700, 'E', 69e9);
mat.steel = struct('name', 'AISI 1080 Steel', 'density', 7850, 'E', 200e9);
mat.carbon_fiber = struct('name', 'Unidirectional Carbon Fiber', 'density', 1600, 'E', 125e9);

mats = mat;
end
