%
% (c) 2016 Charles L. Epstein and Michael O'Neil
%
%     cle@math.upenn.edu
%     oneil@cims.nyu.edu
%
% See the corresponding paper for technical information:
%
%     C. L. Epstein and M. O'Neil, "Smoothed corners and scattered
%     waves", arXiv:1506.08449, 2016.
%     
% This script calls the routines used for smoothing polygons in two
% dimensions. n-gon polygons are simply generated as a list of
% vertices and edges, and then smoothed by convolving with a
% compactly support kernel.
%

addpath('./src');

% specify the number of edges
n0 = 3;

% compute the vertices
ax= [cos(2*pi*(0:(n0-1))/n0),sin(2*pi*(0:n0-1)/n0)];
vx = reshape(ax, [n0 2]);

% set the base rounding size
dh = .05;
n1 = 4;
figure, hold on

for j = 1:n1
    [psa,nt] = smthpoly(vx, j*dh, 4, 2048);
    plot(psa(1:nt,1), psa(1:nt,2), '-', 'LineWidth', 2)
end

% plot the original polygon
nv = sz(vx);
for j =1:nv
    j2 = 1+mod(j,nv);
    X=[vx(j,1),vx(j2,1)];
    Y=[vx(j,2),vx(j2,2)];
    line(X, Y, 'Color', [1,0,0], 'LineWidth', 1);
end

axis equal
hold off
