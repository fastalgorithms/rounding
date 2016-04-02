
addpath('./src');

n0 = 7;
ax= [cos(2*pi*(0:(n0-1))/n0),sin(2*pi*(0:n0-1)/n0)];
vx = reshape(ax, [n0 2]);
dh = .02;
n1 = 8;
figure, hold on
for j = 1:n1
    [psa,nt]=smthpoly(vx, j*dh, 4, 2048);
    plot(psa(1:nt,1),psa(1:nt,2))
end
nv = sz(vx);
for j =1:nv
   
    j2 = 1+mod(j,nv);
    X=[vx(j,1),vx(j2,1)];
     Y=[vx(j,2),vx(j2,2)];
    line(X,Y,'Color',[1,0,0],'LineWidth',2);
end
axis equal
