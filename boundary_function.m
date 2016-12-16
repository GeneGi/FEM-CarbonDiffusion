function [ result ] = boundary_function( x, y, t )

left=-1*pi;
right=1*pi;
bottom=-1*pi;
top=1*pi;

if x==left
    result = 0;
elseif x==right
    result = 0;
elseif y==bottom
    result = -t*sin(x);
elseif y==top
    result = -t*sin(x);
else
    warning='Wrong boundary function or wrong input';
end

end

