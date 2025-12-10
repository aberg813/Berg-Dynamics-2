function [pos, vel] = launcher_ball_state(x, params)
%LAUNCHER_BALL_STATE   Position and velocity of the ball in frame N.
phi = x(1); dphi = x(2); theta = x(3); dtheta = x(4);
LA = params.LA; LB = params.LB;
pos = [LA*cos(phi) + LB*cos(phi + theta);
        LA*sin(phi) + LB*sin(phi + theta)];
vel = [-LA*dphi*sin(phi) - LB*(dphi + dtheta)*sin(phi + theta);
        LA*dphi*cos(phi) + LB*(dphi + dtheta)*cos(phi + theta)];
end
