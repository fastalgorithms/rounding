function [ iout ] = pltsmvx(vp, j0, j1 )
%We plot the smoothed vertex
 nl = sz(vp);
m2 =1;
if nargin == 1
   m0 = 1;
   m1 = nl;
elseif nargin == 2
    m2 = j0;
    m0 = 1;
    m1 = nl;
else
    m0 = j0;
    m1 = j1;
end
for j = m0:m2:m1
    np = sz(vp{j});
    plot3(vp{j}(1:np,1),vp{j}(1:np,2),vp{j}(1:np,3));
end
end

