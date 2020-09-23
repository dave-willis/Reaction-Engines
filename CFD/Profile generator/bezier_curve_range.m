%Range of points in a bezier curve

function bez = bezier_curve_range(n, points)

bez = zeros(n, 2);
for i=1:n
    t = (i-1)/(n-1);
    [bez(i, 1), bez(i, 2)] = bezier(t, points);
end
end