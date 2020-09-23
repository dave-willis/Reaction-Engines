%Calculate coordinate of the point on a bezier curve

function [x, y] = bezier(t, points)

n = length(points) - 1;
x = 0;
y = x;
for i=1:length(points)
    bern = bernstein(t, i-1, n);
    x = x + points(i, 1)*bern;
    y = y + points(i, 2)*bern;
end
end