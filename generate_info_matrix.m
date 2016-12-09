function [P, T, Pb, Tb] = generate_info_matrix(omega, h, basis_type)
% generate_info_matrix
%     - generate the info matrix P, T, Pb, Tb for ordinary mesh, we may
%       just load the info matrix from another software like FreeFEM++
% Input:
%     - (1x4 matrix) omega: [left, right, bottom, top]
%     - (2x1 matrix) h: step size [h_x; h_y]
%     - (string) basis_type: 'linear' or 'quadratic'
% Output:
%     - (2x((Nx+1)*(Ny+1)) matrix) P: the coordinates of all mesh nodes
%     - (3x(Nx * Ny) matrix) T: the global node indices of the mesh nodes
%     - P, T depend on the mesh only, Pb, Tb depend on the type of fem and mesh
%     - Pb: the coordinates of all finite element nodes
%     - Tb: the global node indices of the finite element nodes
%     - for linear basis type, Pb = P, Tb = T
%     - for quadratic basis type, TODO
left = omega(1);
right = omega(2);
bottom = omega(3);
top = omega(4);
hx = h(1);
hy = h(2);
Nx = (right - left) / hx;
Ny = (top - bottom) / hy;
N = (Nx+1) * (Ny+1);

P = zeros(2, N);
T = zeros(3, 2*Nx*Ny);
Q = zeros(Nx+1, Ny+1); %Q is a temp matrix for generating T

for i = 1 : Nx+1
    for j = 1 : Ny+1
        P(1, (i-1)*(Ny+1) + j) = left + (i-1) * hx;
        P(2, (i-1)*(Ny+1) + j) = bottom + (j-1) * hy;
        Q(i, j) = (i - 1) * (Ny + 1) + j;
    end
end

for n = 1 : Nx*Ny
    if mod(n, Ny) == 0
        row = Ny;
        column = n / Ny;
    else
        row = mod(n, Ny);
        column = fix(n / Ny) + 1;
    end

    T(1, 2*n-1) = Q(column, row);
    T(2, 2*n-1) = Q(column+1, row);
    T(3, 2*n-1) = Q(column, row+1);

    T(1, 2*n) = Q(column, row+1);
    T(2, 2*n) = Q(column+1, row);
    T(3, 2*n) = Q(column+1, row+1);
end

if strcmp(basis_type, 'linear')
    Pb = P;
    Tb = T;
elseif strcmp(basis_type, 'quadratic')
    %TODO
else
    error('unknown basis type:' + basis_type)
end

end
