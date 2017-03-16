omega = [-1 1 -1 1];
h = [1/4; 1/4];
basis_type = 'linear';

[P, T, Pb, Tb] = generate_info_matrix(omega, h, basis_type);
[N, Nm, Nb, Nlb] = generate_num(P, T, Pb, Tb);

A1 = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 0);
A2 = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
A3 = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 1, 0);
t_min = 0;
t_max = 3600;
dt = 1800;
tempature = 1183;
[ D ] = func_D(tempature);
A = A1 + dt * D * (A2+A3);
k = 1;
b = assemble_vector_b('func_f', P, T, Pb, Tb, N, Nb, Nlb, t_min+k*dt, D, 0, 0);

vertices = P(:, T(:, 1))
size(Pb)
Nt = (t_max - t_min) / dt;
C = zeros ( Nb, Nt );
