function [ C, maxerror ] = carbon_diffusion_2d( omega, h, t_min, t_max, dt, basis_type, tempature)
% carbon_diffusion_2d:
%     - main function for solving the carbon diffusion equation with FEM method
% Input 
%     - (4x1 matrix) omega: [left, right, bottom, top]
%     - (1x2 matrix) h: step size, [h_x; h_y] 
%     - (scalar) t_min, t_max: the diffusion time
%     - (scalar) dt: time step size
%     - (string) basis_type: 'linear' or 'quadratic'
%     - (scalar) tempature: the diffusion tempature, related to D
% Output
%     - (Nb x Nt matrix) C: the numerical solution for the equation
%     - (scalar) maxerror: compare C with the exact value, get the maximum difference

[P, T, Pb, Tb] = generate_info_matrix(omega, h, basis_type);
[ boundary_nodes ] = generate_boundarynodes(omega, h, basis_type);
[N, Nm, Nb, Nlb] = generate_num(P, T, Pb, Tb);

[ A1 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 0);
[ A2 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
[ A3 ] = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 1, 0);
[ D ] = func_D(tempature);

A = A1 + dt * D * (A2+A3);

[ C ] = time_iteration(P, T, Pb, Tb, boundary_nodes, N, Nb, Nlb, t_min, t_max, dt, A, A1, D)
%[ maxerror ] = get_maxerror('func_exact', C, Pb);
maxerror = 1;
exact = func_exact(P(1, :), P(2, :), 1800)';

end

