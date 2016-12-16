function [ boundary_nodes ] = generate_boundarynodes( omega, h, basis_type )
left = omega(1);
right = omega(2);
bottom = omega(3);
top = omega(4);
hx = h(1);
hy = h(2);
Nx = (right - left) / hx;
Ny = (top - bottom) / hy;
if strcmp(basis_type, 'quadratic')
    Nx = 2 * Nx;
    Ny = 2 * Ny;
end

nbn=2*(Nx+Ny);
boundary_nodes=zeros(2,nbn);

%All Dirichlet boundary nodes.
boundary_nodes(1,:)=-1;

%bottom boundary nodes.
for k=1:Nx
    boundary_nodes(2,k)=(k-1)*(Ny+1)+1;
end

%right boundary nodes.
for k=Nx+1:Nx+Ny
    boundary_nodes(2,k)=Nx*(Ny+1)+k-Nx;
end

%top boundary nodes.
for k=Nx+Ny+1:2*Nx+Ny
    boundary_nodes(2,k)=(2*Nx+Ny+2-k)*(Ny+1);
end

%left boundary nodes.
for k=2*Nx+Ny+1:nbn
    boundary_nodes(2,k)=2*Nx+2*Ny+2-k;
end

end

