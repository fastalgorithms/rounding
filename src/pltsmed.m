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

function [ iout ] = pltsmed(smed,m)
%
% this function plots the smoothed edges of a polyhedron
%

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

