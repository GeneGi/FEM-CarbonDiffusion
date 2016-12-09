function [ D ] = func_D( tempature )
% func_D
%     - caculate the diffusion factor D
%     - D = D0exp(-Q/RT)

D0 = 0.162 * 100;
Q = 137800;
R = 8.314;
T = tempature;
D = D0 * exp( -Q / (R * T));

end

