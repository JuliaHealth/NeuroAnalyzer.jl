function [angle_r, angle_d] = vec_angle(x, y)
	% VEC_ANGLE(x, y)
	% Calculates the angle between two vectors `x` and `y`.

    dotxy = dot(x, y);
    normx = norm(x);
    normy = norm(y);
    angle_r = dotxy / (normx * normy);
    angle_d = (angle_r * 180) / pi;
end