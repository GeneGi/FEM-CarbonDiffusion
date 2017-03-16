omega = [-1 1 -1 1];
h = [1 1];
basis_type = 'linear';

[P, T, Pb, Tb] = generate_info_matrix(omega, h, basis_type);
[boundary_nodes] = generate_boundarynodes(omega, h, basis_type)
[N, Nm, Nb, Nlb] = generate_num(P, T, Pb, Tb);

A = assemble_matrix_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
