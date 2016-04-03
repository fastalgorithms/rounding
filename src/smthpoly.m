%
% (c) 2016 Charles L. Epstein and Michael O'Neil
%
% See gitlab.com/oneilm/rounding/license.md for copyright information,
% and the corresponding paper for technical information:
%
%     C. L. Epstein and M. O'Neil, "Smoothed corners and scattered
%     waves", arXiv:1506.08449, 2016.
%     

function [psamp,ntot] = smthpoly(v,he,ke,me)
%
% This function smooths a polygon
%
% Input variables
%     v are the vertices
%     he, ke are smoothing parameters
%     2*me is the number of sample points are each smoothed vertex
%
% Output variables
%     psamp are sample points lying on the smoothed polygon
%     ntot = sz(psamp)
%
    
    % How many vertices.
    nv = sz(v);

    % We now set up the local geometries for the edges
    for j = 1:nv
        j1 = 1+mod(j,nv);
        j0 =1+mod(j-2,nv);
        X(j,1:2) = (v(j1,1:2)-v(j,1:2));
        X(j,1:2)=X(j,1:2)/sqrt(X(j,1:2)*X(j,1:2)');
        Y(j,1:2) = (v(j0,1:2)-v(j,1:2));
        Y(j,1:2)=Y(j,1:2)/sqrt(Y(j,1:2)*Y(j,1:2)');
    end

    % We get samples of the function used to smooth the vertices

    % We use the part with indices from j0 to j1+1
    [as0,j0,j1]= smthabv(he,ke,2*me);
    as = as0(j0:j1+1);
    s0 = (0:4*me-1)/(2*me)-1;
    s = s0(j0:j1+1);
    x = as+s;
    y = as-s;
    
    % The number of points on the smoothed edges.
    npte = j1-j0+2; 

    % Find the intersections of the smoothed edges near to the vertex
    % A tolerance for calculations

    % The number of sample points on all the curved parts.
    npt1 = nv*npte; 
    nx2 = ceil(log2(npt1));
    nflp = nx2;

    % Loop over vertices
    for q2 = 1:nv 

        psamp(1+(q2-1)*(npte+nflp):q2*npte+(q2-1)*nflp,1:2) = ...
            ones(npte,1)*v(q2,1:2)+(x(1:npte))'*X(q2,1:2) ...
            +(y(1:npte))'*Y(q2,1:2);

        % The last point on the current edge  
        Z1 = v(q2,1:2)+x(npte)*X(q2,1:2)+y(npte)*Y(q2,1:2); 

        % Next q2 index in cyclic order
        q2n=1+mod(q2,nv); 

        % The first point on the next edge
        Z2 = v(q2n,1:2)+x(1)*X(q2n,1:2)+y(1)*Y(q2n,1:2);
            
        % Test code:
        Zdif = Z2-Z1;
        nZ =sqrt(Zdif*Zdif');

        % Construct the straight segment joining the two curved ones
        psamp(q2*npte+(q2-1)*nflp+1:q2*(npte+nflp),1:2) = ...
            ones(nflp,1)*Z1(1:2)+((1:nflp)'/(1+nflp))*Zdif(1:2);
    end
    ntot = sz(psamp);
    psamp(ntot+1,1:2)=psamp(1,1:2);
    ntot = ntot+1;

    % To get a closed curve
    %psamp(ntot+1)=psamp(1);
end
  

