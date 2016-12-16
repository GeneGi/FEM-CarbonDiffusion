function [ A ] = assemble_matrix_A( P, T, Pb, Tb, N, Nb, Nlb, derative_degree_x, derative_degree_y)
% assemble_matrix_A:
%     - 
% Input:
%     - (matrix) Pb: the coordinates of all finite element nodes
%     - (matrix) Tb: the global node indices of the finite element nodes
%     - (scalar) N, Nb, Nlb: number of mesh elements, mesh nodes and local fem nodes
%     - (scalar) derative_degree_x, derative_degree_y: derative degree of basis function
% Output:
%     - (Nb x Nb matrix) A: the stiffness matrix
A = sparse(Nb,Nb);

for n = 1 : N
    vertices = P(:, T(:, n));
    for alpha = 1 : Nlb
        for beta = 1 : Nlb
            temp = gauss2d_integral_trial_test(vertices, alpha, beta, derative_degree_x, derative_degree_y);
            A(Tb(beta, n),Tb(alpha, n)) = A(Tb(beta,n),Tb(alpha,n)) + temp;
        end
    end
end

end

