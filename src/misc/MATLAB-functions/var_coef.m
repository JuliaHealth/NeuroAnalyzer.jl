function var_coef = var_coef(se, p)
	% VAR_COEF(se, p)
	% Coeffcient of variation.
	% Calculates coefficient of variation for standard error se and statistics p (e.g. mean).

    var_coef =  100 * (se / p);
end