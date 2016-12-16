function [ A, b ] = treat_boundary( boundary_func, boundary_nodes, A, b, Pb, t)

nbn = size(boundary_nodes,2);

for k = 1:nbn
    if boundary_nodes(1,k) == -1
        i = boundary_nodes(2,k);
        A(i,:) = 0;
        A(i,i) = 1;
        b(i) = feval(boundary_func,Pb(1,i),Pb(2,i), t);
    end
end


end

