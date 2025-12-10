function [value, isterminal, direction] = launcher_release_event(~, x, params)
%LAUNCHER_RELEASE_EVENT  Event for ode45: x_ball = 0 indicates release.
[pos, ~] = launcher_ball_state(x, params);
value = pos(1);
isterminal = 1;
direction = -1;
end
