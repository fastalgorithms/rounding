function [isin] = asubb(A,B)
%
% This function compares two lists of identical types of data
% (e.g. integers) and returns 1 if A is a subset of B, and zero
% otherwise.  It does not matter if elements of A or B are repeated,
% they are effectively treated as sets.
%
    na= max(size(A));
    nb= max(size(B));
    ncnt=0;
    for j = 1:na
        ninc = 0;
        for k = 1:nb
            if B(k) == A(j) ninc =1;
            end
        end
        ncnt = ncnt + ninc;
    end
    if ncnt == na 
        isin = 1;
    else
        isin = 0;
    end
end

