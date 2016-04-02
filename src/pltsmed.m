function [ iout ] = pltsmed(smed,m)
%We plot the smoothed edge
nl = sz(smed);
m0 = 1;
if nargin == 2
    m0 =m;
end
for j = 1:m:nl
    np = sz(smed{j});
    plot3(smed{j}(1:np,1),smed{j}(1:np,2),smed{j}(1:np,3));
end
iout = 1;
end

