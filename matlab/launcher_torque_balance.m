function val = launcher_torque_balance(theta, params)
%LAUNCHER_TORQUE_BALANCE  Static torque balance for theta when phi = 0.
val = params.LB * params.g * (0.5*params.mB + params.mBall) * cos(theta) + params.k * theta;
end
