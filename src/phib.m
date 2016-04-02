function  [x]= phib( b,t )
%
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
    x=t.^(1-b).*besselj(b-1,2*sqrt(1+b)*t);
end

