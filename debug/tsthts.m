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
% merely some testing code for the polygonal convolutions
%

addpath('../src')

overlap = 5;
m=512;
hv = .05;
kv = 12;

% We use the part with indices from 2m+1 to 2m+1+(j1-j0+1)/2 + overlap
[as1,j2,j3]= smthabv(hv,kv,2*m);
s1 = (0:4*m-1)/(2*m)-1; 
mlst = 2*m+1+((j3-j2+1)/2)+overlap;

% This function starts linear and ends up smoothed
bs = as1(mlst:-1:2*m+1); 

% This function starts positive and ends up at 0.
ss = abs(s1(mlst:-1:2*m+1)); 

figure, plot(s1,as1)
figure, plot(ss,bs)

% This is the number of heights we will use (this could be decimated later)
nhts = sz(ss) 
