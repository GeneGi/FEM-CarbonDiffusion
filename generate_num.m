function [N, Nm, Nb, Nlb] = generate_num(P, T, Pb, Tb)
% generate_num:
%     - generate the basic num for info matrix
% Input:
%     - (2x((Nx+1)*(Ny+1)) matrix) P: the coordinates of all mesh nodes
%     - (3x(Nx * Ny) matrix) T: the global node indices of the mesh nodes
%     - P, T depend on the mesh only, Pb, Tb depend on the type of fem and mesh
%     - Pb: the coordinates of all finite element nodes
%     - Tb: the global node indices of the finite element nodes
% Output:
%     - (scalar) N: number of mesh elements
%     - (scalar) Nm: number of mesh nodes
%     - (scalar) Nb: number of finite element nodes
%     - (scalar) Nlb: number of local finite element nodes in a mesh element

N = size(T, 2);
Nm = size(P, 2);
Nlb = size(Tb, 1);
Nb = size(Pb, 2);

end

