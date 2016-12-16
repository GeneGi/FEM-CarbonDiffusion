function [ C, E ] = carbon_diffusion_2d( omega, h, Nt,  dt, basis_type, tempature)
% carbon_diffusion_2d:
%     - main function for solving the carbon diffusion equation with FEM method
% Input 
%     - (4x1 matrix) omega: [left, right, bottom, top]
%     - (1x2 matrix) h: step size, [h_x; h_y]
%     - (scalar) Nt: number of time
%     - (scalar) dt: time step size
%     - (string) basis_type: 'linear' or 'quadratic'
%     - (scalar) tempature: the diffusion tempature, related to D
% Output
%     - (Nb x Nt matrix) C: the numerical solution for the equation
%     - (Nb x Nt matrix) E: compare C with the exact value, get the maximum difference

[ P, T, Pb, Tb ] = generate_info_matrix(omega, h, basis_type);
[ boundary_nodes ] = generate_boundarynodes(omega, h, basis_type);
[ N, Nm, Nb, Nlb ] = generate_num(P, T, Pb, Tb);

[ A1 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 0);
[ A2 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
[ A3 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 1, 0);
[ D ] = func_D(tempature);

A = A1 + dt * D * (A2+A3);

C = zeros ( Nb, Nt );
C(:, 1) = initial_function(Pb(1, :), Pb(2, :));
for k = 1 : Nt-1
    b1 = assemble_vector_b('func_f', P, T, Pb, Tb, N, Nb, Nlb, k*dt, D, 0, 0);
    b2 = A1 * C(:, k);
    b = dt * b1 + b2;
    [A, b] = treat_boundary('boundary_function', boundary_nodes, A, b, Pb, dt*k);
    C(:, k+1) = A \ b;
end

for k = 1 : Nt
    E(:, k) = get_maxerror('func_exact', C(:, k), Pb, dt*(k-1));
end

