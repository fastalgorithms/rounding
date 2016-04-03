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

function [samps,j0,j1] = smthabv(h,k,n)
%
% Using the FFT this function returns 2n samples of the convolution
% between |x| and (1-(x/h)^2)^k*chi_{[-1,1]}(x). Note the minimum
% value of the smoothed function occurs at index n+1 We also return
% the indices of the endpoints of the curved section, which coincides
% with the intersections of the smoothed curve with the faces.
%
    
%    
% We first compute the Fourier coefficients of |x| exactly
%

    % zero frequency
    ah(1)=pi^2; 
    
    % the positive frequencies
    ah(2:n) = -2*(1-(-1).^(1:(n-1)))./(1:(n-1)).^2; 
    
    % the negative frequencies
    ah(n+1:2*n) = -2*(1-(-1).^(n:-1:1))./(n:-1:1).^2; 

    % We comnpute samples of the (unscaled) smoothing kernel
    j0 = floor(n*(1-h)+1);
    j1 = floor(n*(1+h)+1);

    % We take account of the characteristic function of the interval
    % [-h,h]:
    fk(1:j0)=0;
    fk(j0+1:j1)=(1-((j0:(j1-1))./(n*h)-1/h).^2).^k;
    fk(j1+1:2*n)=0;
    ftfk = fft(fk);
    
    % Multiply the ft of |x| by the ft of the filter
    filftfk = ah.*ftfk; 
    
    % Invert to get the spatial samples.
    samps = n*real(ifft(filftfk))/(pi^2*ftfk(1));

end
