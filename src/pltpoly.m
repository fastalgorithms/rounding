%
% (c) 2016 Charles L. Epstein and Michael O'Neil
%
% See gitlab.com/oneilm/rounding/license.md for copyright information,
% and the corresponding paper for technical information:
%
%     C. L. Epstein and M. O'Neil, "Smoothed corners and scattered
%     waves", arXiv:1506.08449, 2016.
%     

function [iout] = pltpoly(V,E)
%
% This function plots a polyhedron, assuming that the figure has
% been created outside of the function
%
%   V is a list of vertices [v1;v2;...;vn]
%   E is a set of edges {[i1,i2],...[j1,j2]}
%
    ne = sz(E);
    X=[V(E{1}(1),1),V(E{1}(2),1)];
    Y=[V(E{1}(1),2),V(E{1}(2),2)];
    Z=[V(E{1}(1),3),V(E{1}(2),3)];

    line(X,Y,Z,'Color',[1,0,0],'LineWidth',1);
    for j=2:ne
        X=[V(E{j}(1),1),V(E{j}(2),1)];
        Y=[V(E{j}(1),2),V(E{j}(2),2)];
        Z=[V(E{j}(1),3),V(E{j}(2),3)];
        line(X,Y,Z,'Color',[1,0,0],'LineWidth',1);
    end
    iout = 1;
end

