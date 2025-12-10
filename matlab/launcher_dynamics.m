function dx = launcher_dynamics(~, x, params)
%LAUNCHER_DYNAMICS   State derivative for the two-link launcher.
%   State vector x = [phi, phidot, theta, thetadot].

phi = x(1); dphi = x(2); theta = x(3); dtheta = x(4);
[M, forcing] = launcher_mass_forcing(phi, theta, dphi, dtheta, params);
acc = M \ forcing;

dx = zeros(4,1);
dx(1) = dphi;
dx(2) = acc(1);
dx(3) = dtheta;
dx(4) = acc(2);
end

function [M, forcing] = launcher_mass_forcing(phi, theta, dphi, dtheta, params)
LA = params.LA; LB = params.LB;
mA = params.mA; mB = params.mB; mBall = params.mBall;
IA = params.IA; IB = params.IB;
k = params.k; g = params.g; tau = params.tau;

cth = cos(theta);
sth = sin(theta);
cphi = cos(phi);
cphitheta = cos(phi + theta);

i11 = IA + IB + 0.25*LA^2*mA + 0.25*mB*(4*LA^2 + 4*LA*LB*cth + LB^2) + ...
       mBall*(LA^2 + 2*LA*LB*cth + LB^2);
i12 = IB + 0.25*LB*mB*(2*LA*cth + LB) + LB*mBall*(LA*cth + LB);
i22 = IB + 0.25*LB^2*mB + LB^2*mBall;

M = [i11, i12; i12, i22];

forcing1 = 0.5*LA*LB*mB*dtheta*(2*dphi + dtheta)*sth + ...
           LA*LB*mBall*dtheta*(2*dphi + dtheta)*sth - 0.5*LA*g*mA*cphi - ...
           0.5*g*mB*(2*LA*cphi + LB*cphitheta) - g*mBall*(LA*cphi + LB*cphitheta) + tau;
forcing2 = -0.5*LA*LB*mB*dphi^2*sth - LA*LB*mBall*dphi^2*sth - ...
           0.5*LB*g*mB*cphitheta - LB*g*mBall*cphitheta - k*theta;

forcing = [forcing1; forcing2];
end
