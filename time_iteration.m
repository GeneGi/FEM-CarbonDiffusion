function [ C ] = time_iteration(P, T, Pb, Tb, boundary_nodes, N, Nb, Nlb, t_min, t_max, dt, A, A1, D)

Nt = (t_max - t_min) / dt;
C = zeros ( Nb, Nt );
C(:, 1) = initial_function(Pb(1, :), Pb(2, :));
for k = 1 : Nt-1
    b1 = assemble_vector_b('func_f', P, T, Pb, Tb, N, Nb, Nlb, t_min+k*dt, D, 0, 0);
    b2 = A1 * C(:, k);
    b = dt * b1 + b2;
    [A, b] = treat_boundary('boundary_function', boundary_nodes, A, b, Pb, dt*(k-1));
    C(:, k+1) = A \ b;
end

end
