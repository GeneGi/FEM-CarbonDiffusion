function [ quad ] = gauss2d_integral_trial_test(vertices, basis_index_trial, basis_index_test, derative_degree_x, derative_degree_y)
x1 = vertices(1, 1);
y1 = vertices(2, 1);
x2 = vertices(1, 2);
y2 = vertices(2, 2);
x3 = vertices(1, 3);
y3 = vertices(2, 3);

area = abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;

A = [ 0.33333333333333 0.47014206410511 0.47014206410511 0.05971587178977 0.10128650732346 0.10128650732346 0.79742698535309;
      0.33333333333333 0.47014206410511 0.05971587178977 0.47014206410511 0.10128650732346 0.79742698535309 0.10128650732346 ];
W = [ 0.22500000000000 0.13239415278851 0.13239415278851 0.13239415278851 0.12593918054483 0.12593918054483 0.12593918054483 ];

quad = 0.0;
Ng = length(W);

for i = 1 : Ng
    x = x1 * ( 1 - A(1, i) - A(2, i)) + x2 * A(1, i) + x3 * A(2, i);
    y = y1 * ( 1 - A(1, i) - A(2, i)) + y2 * A(1, i) + y3 * A(2, i);
    quad = quad + W(i) * basis_function(x, y, vertices, basis_index_trial, derative_degree_x, derative_degree_y)...
                       * basis_function(x, y, vertices, basis_index_test, derative_degree_x, derative_degree_y);
end
quad = area * quad;

end

