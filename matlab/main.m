function main
% Finite Element Method for solving 2d Carbon Diffusion Equation
%      Ct = D \nabla C
% C is the Carbon Density
% D is the Diffusion factor which is constant here

H = [1/4 1/8 1/16 1/32 ;    %step size x
     1/4 1/8 1/16 1/32 ];   %step size y

H_num = size(H, 2);

E = zeros(H_num, 1);        %max error
Order = zeros(H_num-1, 1);

omega = [-1 1 -1 1]; %[left, right, bottom, top]

t_min = 0;
t_max = 3600;
dt = 1800;               %time step size

basis_type = 'linear';
tempature = 1183;

for i = 1 : H_num
    h = H(:, i);
    [C, E(i)] = carbon_diffusion_2d(omega, h, t_min, t_max, dt, basis_type, tempature);
end

for i = 1 : (H_num-1)
    Order(i) = log(E(i)/E(i+1))/log(H(1,i)/H(1,i+1));
end

end
