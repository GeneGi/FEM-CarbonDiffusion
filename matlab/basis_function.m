function [ result ] = basis_function( x, y, vertices, basis_index, derative_degree_x, derative_degree_y )
x1 = vertices(1, 1);
y1 = vertices(2, 1);
x2 = vertices(1, 2);
y2 = vertices(2, 2);
x3 = vertices(1, 3);
y3 = vertices(2, 3);
J = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

if derative_degree_x == 0 && derative_degree_y == 0
    if basis_index == 1
        result = - ((y3 - y1) * (x - x1) - (x3 - x1) * (y - y1)) / J...
                 - (-(y2 - y1) * (x - x1) + (x2 - x1) * (y - y1)) / J + 1;
    elseif basis_index == 2
        result = ((y3 - y1) * (x - x1) - (x3 - x1) * (y - y1)) / J;
    elseif basis_index == 3
        result = (-(y2 - y1) * (x - x1) + (x2 - x1) * (y - y1)) / J;
    else
        error('Error: unknown basis index.\n');
    end
elseif derative_degree_x == 1 && derative_degree_y == 0
    if basis_index == 1
        result = (y2 - y3) / J;
    elseif basis_index == 2
        result = (y3 - y1) / J;
    elseif basis_index == 3
        result = -(y2 - y1) / J;
    else
        error('Error: unknown basis index.\n');
    end
elseif derative_degree_x == 0 && derative_degree_y == 1
    if basis_index == 1
        result = (x3 - x2) / J;
    elseif basis_index == 2
        result = -(x3 - x1) / J;
    elseif basis_index == 3
        result = (x2 - x1) / J;
    else
        error('Error: unknown basis index.\n');
    end
else
    error('Error: unknown derative degree.\n');
end

