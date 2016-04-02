function [edgesamp] = smthedge(v0,v1,nv0,nv1,mv0,mv1,lmin0,lmin1,h,k,ne,nh)
%This function smooths an edge which joins v0 to v1 where
%nv0 and nv1 are the outer normal vectors to the faces that interesect
%along the given edge. h and k are the smoothing parameters.
%lmin0 and lmin1 are the "top" and "bottom" of the smoothed vertex in
%terms of <X-v0,mv0> and <X-v1,mv1> respectively
%The program outputs 2*n+1 samples of the smoothed edge 
%lying in the plane orthogonal to the edge passing through its midpoint.
%   
m = (v0+v1)/2; %the midpoint of the edge
w = (v0-v1)/(sqrt((v0-v1)*(v0-v1)'));% the unit direction from v0 to v1
if det([w;nv0;nv1])>0 %We find vectors normal to the edge lying along the two 
                      %faces that meet along the edge.
    X = cp(nv0,w);
    Y = cp(w,nv1);
else
    X = cp(nv1,w);
    Y = cp(w,nv0);
end
[as0,j0,j1] = smthabv(h,k,2*ne);
s0 = (0:4*ne-1)/(2*ne)-1;
as = as0(j0:j1);
s = s0(j0:j1);
%Projections into the top and bottom planes
xp0 = X*mv0';
yp0 = Y*mv0';
xp1 = X*mv1';
yp1 = Y*mv1';
wp0 = w*mv0';
wp1 = w*mv1';
nj = j1-j0+1;
%Now we produce the samples along the smoothed edge
for j2 = 1:nj
     x = (as(j2)+s(j2));
     y = (as(j2)-s(j2));
    tmax = (lmin0-(m-v0)*mv0'-x*xp0-y*yp0)/wp0;
    tmin = (lmin1-(m-v1)*mv1'-x*xp1-y*yp1)/wp1;
    tt = tmax*(1-((0:(nh-1))/(nh-1)))+tmin*(0:(nh-1))/(nh-1);
edgesamp{j2}(1:nh,1:3) = ones(nh,1)*(m(1:3)+x*X(1:3)+y*Y(1:3))+tt(1:nh)'*w(1:3);
end
end

