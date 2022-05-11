T = zeros(3,3);
l = zeros(3,1);
for i = 1:iter
    T = T + [y(1,i); y(2,i); y(3,i)]*[y(1,i); y(2,i); y(3,i)].'
    l = l+ [y(1,i); y(2,i); y(3,i)]
end;

rank(T-1/iter*l*l.')

B = [T l;l.' iter];

rank(B)

det(B)

l