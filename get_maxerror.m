function [ maxerror ] = get_maxerror(exact_function, C, Pb, t)

error = C - feval(exact_function, Pb(1, :), Pb(2, :), t)';
maxerror = max(abs(error));

end
