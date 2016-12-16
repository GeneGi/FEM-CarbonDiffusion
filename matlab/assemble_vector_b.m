function [ b ] = assemble_vector_b( func_name, P, T, Pb, Tb, N, Nb, Nlb, time, D, derative_degree_x, derative_degree_y) 
b = zeros(Nb,1);
for n = 1 : N
    vertices = P(:,T(:,n));
    for beta = 1 : Nlb
        temp = gauss2d_integral_test(func_name, vertices, beta, time, D, derative_degree_x, derative_degree_y);
        b(Tb(beta,n),1) = b(Tb(beta,n),1) + temp;
    end
end

end

