%
% (c) 2016 Charles L. Epstein and Michael O'Neil
%
% See gitlab.com/oneilm/rounding/license.md for copyright information,
% and the corresponding paper for technical information:
%
%     C. L. Epstein and M. O'Neil, "Smoothed corners and scattered
%     waves", arXiv:1506.08449, 2016.
%     
% This script calls the routines used for smoothing polyhedra. It
% begins with data tables describing standard platonic solids embedded
% in 3-space as a list of vertices, V; and then the indices of
% vertices defining the edges, E; and lists of vertices defining the
% faces, F.
%

addpath('./src');


%
% The following built-in polyhedrons are available, and determined
% by the value of pnum:
%
%   1 - pyramid
%   2 - cube
%   3 - regular tetrahedron
%   4 - octohedron
%   5 - icosahedron
%

%
% First choose a polyhedron
%
pnum = 5


%
% Load the requistie vertices and edges
%

if pnum == 1
    % pyramid
    V = [[0,0,1];[-1,0,0];[0,1,0];[1,0,0];[0,-1,0]]
    E = {[1,2],[1,3],[1,4],[1,5],[2,3],[3,4],[4,5],[5,2]}
    F = {[1,2,3],[1,3,4],[1,4,5],[1,5,2],[2,3,4,5]}
elseif pnum == 2
    % Data for a cube
    V = [[0,0,0];[1,0,0];[1,1,0];[0,1,0];[0,1,1];[0,0,1];...
         [1,0,1];[1,1,1]]
    E = {[1,2],[2,3],[3,4],[4,1],[1,6],[2,7],[3,8],[4,5],[5,6],...
         [6,7],[7,8],[8,5]}
    F = {[1,2,3,4],[1,2,7,6],[1,4,5,6],[3,4,5,8],[2,3,7,8],[5,6,7,8]}
elseif pnum ==3
    % Data for a regular tetrahedron
    V = [[1,0,0]; [cos(2*pi/3),sin(2*pi/3),0]; [cos(4*pi/3),...
             sin(4*pi/3),0];[0,0,sqrt(2-2*cos(2*pi/3))]]
    E = {[1,2],[2,3],[3,1],[4,1],[4,2],[4,3]}
    F = {[1,2,3],[1,2,4],[2,3,4],[1,3,4]}
elseif pnum == 4
    % Data for an octohedron
    V = [[0,0,1];[-1,0,0];[0,1,0];[1,0,0];[0,-1,0];[0,0,-1]]
    E = {[1,2],[1,3],[1,4],[1,5],[2,3],[3,4],[4,5],[5,2],...
         [2,6],[3,6],[4,6],[5,6]}
    F = {[1,2,3],[1,3,4],[1,4,5],[1,5,2],[3,4,6],[4,5,6],...
         [5,2,6],[2,3,6]}
elseif pnum == 5
    % Data for an icosahedron
    V= [[0, (1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4), ...
         (1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5))];...
        [0, (1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4),...
         -(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5))];...
        [0, (1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5)),...
         (1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5))];...
        [0, (1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5)),...
         -(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5))];...
        [(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0, ...
         (1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4)];...
        [(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0, ...
         (1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5))];...
        [-(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0, ...
         (1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4)];...
        [-(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0, ...
         (1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5))];...
        [(1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4), ...
         (1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0];...
        [(1/5)*sqrt(1/2+(1/2)*sqrt(5))*5^(3/4), ...
         -(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0];...
        [(1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5)), ...
         (1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0];...
        [(1/5)*5^(3/4)*(-1/2-(1/2)*sqrt(5))/sqrt(1/2+(1/2)*sqrt(5)), ...
         -(1/5)*5^(3/4)/sqrt(1/2+(1/2)*sqrt(5)), 0]]
    E= {[7,5],[7,1],[7,11],[7,12],[7,3],[5,1],[5,9],[5,10],[5,3],...
        [9,1],[9,2],[9,10],[9,6],[10,4],[10,6],[10,3],[6,2],[6,8],[6,4],...
        [4,8],[4,12],[4,3],[3,12],[1,2],[1,11],[2,8],[2,11],[11,12],...
        [11,8],[8,12]}
    F= {[2,6,9],[2,9,1],[9,1,5],[2,1,11],[2,11,8],...
        [8,11,12],[8,12,4],[8,4,6],...
        [4,12,3],[4,3,10],[4,10,6],[6,10,9],[9,10,5],[10,3,5],...
        [5,3,7],[3,7,12],...
        [7,11,12],[1,7,11],[8,6,2],[7,1,5]}
end

%
% Obtain data structures used in the smoothing procedure
%
[V0,E0,F0,lm] = polyframe(V,E,F);


%
% lm is the length of the shortest edge. It establishes a scale for all
% subsequent calculations. (the smoothing parameter "he" can range from 
% 0.0000001*lm up to about 0.3*lm.)
%
he = .005*lm; 
ke = 6; 
hv =6*he*lm; 
kv = 4; 
me=512; 
mv=512; 
nh=512;

figure
axis equal
axis off
hold on

%
% We smooth the vertices and plot the result
%
for j =1:sz(V)
    [smvp, lmin(j)] = smthvrtx(j,V,E,V0,F0,he,ke,hv,kv,me,mv);
    pltsmvx(smvp,10);
end

%
% We smooth the edges and plot the result
%
for j = 1:sz(E)

    % Get the indices of the ends of the edge
    i=E0{j,6};
    v0 = V(i(1),1:3);
    v1 = V(i(2),1:3);

    % The normal vectors of the planes defining the edge
    nv0 = E0{j,4}(1:3);
    nv1 = E0{j,5}(1:3);

    % Get the outer support vectors at the vertices
    mv0 = V0{i(1),2};
    mv1 = V0{i(2),2};

    % Get the correct lmin data
    lmin0 = lmin(i(1));
    lmin1 = lmin(i(2));

    % Get the points on the smoothed edge
    [edgesamp] = smthedge(v0,v1,nv0,nv1,mv0,mv1,lmin0,lmin1,he,ke,me,nh);
   pltsmed(edgesamp,2);
end

pltpoly(V,E);