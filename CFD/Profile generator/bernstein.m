%Bernstein polynomial

function bern = bernstein(t, i, n)

bern = nchoosek(n, i)*(t^i)*(1-t)^(n-i);

end


