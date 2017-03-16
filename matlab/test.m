function test
% Finite Element Method for solving 2d Carbon Diffusion Equation
%      Ct = D \nabla C
% C is the Carbon Concentration
% D is the Diffusion factor which is constant here

H = [1/4 1/8 1/16 1/32 ;    %step size x
     1/4 1/8 1/16 1/32 ];   %step size y
H = H * pi;
H_num = size(H, 2);

omega = [-pi, pi, -pi, pi]; %[left, right, bottom, top]

t_min = 0;
t_max = 3600;
dt = 1800;               %time step size
Nt = (t_max - t_min) / dt + 1;

E = zeros(H_num, Nt);        %max error
Order = zeros(H_num-1, 1);  %order of convergence

basis_type = 'linear';
tempature = 1183;        %the Diffusion factor D is related to tempature

%for each step size h, caculate the numerical solution and maxerror
for i = 1 : 3%H_num
    h = H(:, i);
    [C, E(i, :)] = carbon_diffusion_2d(omega, h, Nt, dt, basis_type, tempature)
    
end

%caculate the order of convergence with the max error for each h
for k = 1 : Nt
    for i = 1 : (H_num-1)
        Order(i, k) = log(E(i, k)/E(i+1, k))/log(H(1,i)/H(1,i+1));
    end
end

Order
end
