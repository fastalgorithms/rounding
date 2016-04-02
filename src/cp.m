function [uxv] = cp(u,v)
%
% This function computes the cross product of two 3-vectors
%
    uxv = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3),...
           u(1)*v(2)-u(2)*v(1)];
end

