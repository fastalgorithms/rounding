function [vrtxsamp,lmin] = smthvrtx(i,V,E,Vout,Fout,he,ke,hv,kv,me,mv)
%This function smooths the ith vertex. It does this by using the smoothing
%of the edges to construct smoothed intersections with <X-v0,nv0>=a for
%various values of a. These are extended harmonically and used to smooth
%off the vertex.
%
%Input parameters:
%
%(V,E) a description of the vertices and edges of the original polyhedron.%
%(Vout,Fout) the output of polyframe which provides local descriptions of
%               the geometry in a neighborhood of each vertex and face.
%
%he, ke are smoothing parameters for the edges. 
%
%me is used to construct smoothing function for the edges;
%       it's the number of sample points.
%
%hv, kv are smoothing parameters for the vertex itself. 
%
%mv is used to construct the smoothing function for vertex.
% 
%Outputs:
%
%vrtxsamp = samples of the smoothed vertex.
%lmin = a measure of the size of the smoothed vertex, which is needed to
%       assess how much of the edge to smooth.
% 
%
v0 = V(i,1:3);      % this is the vertex we will smooth.
nv0 = Vout{i,2};    %The unit support vector at the vertex.
ne = sz(Vout{i,1}); %The valence of the vertex.
vind = Vout{i,1};   %The indices of the relevant vertices.
find = Vout{i,3};   %The indices of the corresponding faces.
%
%We now set up the local geometries for the edges that end at the vertex
%we are smoothing.
for j = 1:ne
   v1(1:3) = V(vind(j),1:3);            %Current "other vertex"
   m(j,1:3) = (v0+v1)/2;                %the midpoint of the jth edge
   w(j,1:3) = (v0-v1)/(sqrt((v0-v1)*(v0-v1)'));% the unit direction from v0 to v1
   nv1 = Fout{find(j),1};               %The normal vectors to the faces
   nv2 = Fout{find(1+mod(j,ne)),1};     %defining the current edge.
   if det([w(j,1:3);nv1;nv2])>0 
     X(j,1:3) = cp(nv1,w(j,1:3));       %Local frame along the jth edge
     Y(j,1:3) = cp(w(j,1:3),nv2);
   else
   Y(j,1:3) = cp(nv1,w(j,1:3));         %Local frame along the jth edge
   X(j,1:3) = cp(w(j,1:3),nv2);
   end
end
%We get samples of the function used to smooth the edges
[as0,j0,j1]= smthabv(he,ke,2*me);% We use the part with indices from j0 to j1+1
as = as0(j0:j1+1);
s0 = (0:4*me-1)/(2*me)-1;
s = s0(j0:j1+1);
x = as-s;
y = as+s;
npte = j1-j0+2;             % The number of points on the smoothed edges.
%Find the intersections of the smoothed edges near to the vertex
%A tolerance for calculations
tol = abs(y(1)*100);
for j=1:ne
    j2 = 1+mod(j,ne); % j+1 mod ne
    M2 = (m(j2,1:3)+x(1)*X(j2,1:3));
    M1 = (m(j,1:3)+y(npte)*Y(j,1:3));
    M = M2-M1;
    a = M*w(j,1:3)';
    b = M*w(j2,1:3)';
    c = w(j,1:3)*w(j2,1:3)';
    t1 = (a-b*c)/(1-c*c);
    t2 = (a*c-b)/(1-c*c);
    l(j) = (M2+t2*w(j2,1:3)-v0)*nv0';%The height of the intersection rel nv0.
    if abs(l(j)-(M1+t1*w(j,1:3)-v0)*nv0')> tol 
        error('Problem with intersection of edges %d and %d',j,j2)
    end
    ml(j)=(m(j,1:3)-v0)*nv0'; % The heights of the midpoints.
end
lmax = min(l(1:ne));            % The maximum intersection height rel 
                                % nv0 for smoothing  the vertex.
mmax =  max(ml(1:ne));          % The height of the highest midpoint.                    
%
%We find the 't' coordinates for projection of the smoothed edges into planes
%of the form <X-v0,nv0>=l.
for j=1:ne
    dd(j) = w(j,1:3)*nv0';              %The denominator;
    cc =-((m(j,1:3)-v0)*nv0')/dd(j);    %The constant term;
    xp = -(X(j,1:3)*nv0')/dd(j);
    yp = -(Y(j,1:3)*nv0')/dd(j);
    t(j,1:npte)= (x(1:npte)'*xp+y(1:npte)'*yp)+cc;
end
%We now construct the function for smoothing the vertex in the radial
%direction
overlap = 3;
[ar1,j2,j3]= smthabv(2*hv,kv,2*mv);% We use the part with indices from 
r1 = (0:4*mv-1)/(2*mv)-1;          % 2m+1 to 2m+1+(j1-j0+1)/2 + overlap
mlst = 2*mv+1+((j3-j2+1)/2)+overlap;
ars = ar1(mlst:-1:2*mv+1); %This function starts linear and ends up smoothed
rs = abs(r1(mlst:-1:2*mv+1)); %This function starts positive and ends up at 0.
nhts = sz(rs); % This is the number of heights we will use (this could be decimated later)
%
%Now we construct a function to smooth in the height direction
% For each desired height we now construct a smoothed link of the vertex
[ah1,jh2,jh3]= smthabv(hv,kv,2*mv);% We use the part with indices from 
mdif = nhts-(jh3-jh2+1);
if mdif<= 1
    error('Cannot smooth vertex height, mdif= %d',mdif)
end
hh1 = (0:4*mv-1)/(2*mv)-1; % 
mdif1 = floor(mdif/2);
mdif2 = mdif-mdif1;
%We select points from flat part to flat part.
ahs = ar1(jh2-mdif1:jh3+mdif2); %This function starts linear and ends up smoothed
hs = hh1(jh2-mdif1:jh3+mdif2);
%We fix a spacing for the separation of the intersection planes
nskip = 2; %To compress the smoothed vertex.
dh = abs(mmax)/(nskip*mv);
lmin = lmax-dh*(nhts-1);
ah2 = ahs-hs;
%We need to rescale to get the correct first and last heights
ht(1:nhts) = ((lmin-lmax)*ah2(1:nhts)+...
    ones(1,nhts)*(ah2(1)*lmax-ah2(nhts)*lmin))/(ah2(1)-ah2(nhts));
%Or  use these linear spaced heights...
%We need to change this variable to use a
%collection of heights ht(a) a in [0,1] so the ht(a) is very flat at a=1
%ht(1:nhts)=lmax-dh*(nhts:-1:1);
%We need to assemble the intersection of the planes <X-v0,nv0>=ht as
%closed curves lying in a plane, and then apply the harmonic extension
%technique to smooth vertex.
npt1 = ne*npte; % The number of sample points on all the curved parts.
nx2 = ceil(log2(npt1));
npt2 = 2^(nx2+2)-npt1; % The maximum number of points on all the straight parts;
for q1 =1:nhts %Loop over the set of heights
    nflp = ceil((npt2/ne)*(1-(q1-1)/nhts)); % We reduce the number of points/face on 
          %the flat parts as we get  closer to the vertex. We use the same
          %number of points on each flat section.
          % 
          clear Z Zf Zh
    for q2 = 1:ne %Loop over the edges meeting at the vertex
        %We construct a representation of the intersection of the
        %edge-smoothed polyhedron with the planes <X-v0,nv0>=ht:
        Z(1+(q2-1)*(npte+nflp):q2*npte+(q2-1)*nflp,1:3) = ...
            ones(npte,1)*m(q2,1:3)+(x(1:npte))'*X(q2,1:3)+ ...
            (y(1:npte))'*Y(q2,1:3)+...
            (t(q2,1:npte)'+(ht(q1)/dd(q2))*ones(npte,1))*w(q2,1:3);
        Z1 = m(q2,1:3)+x(npte)*X(q2,1:3)+ ...%The last point on the current edge
            y(npte)*Y(q2,1:3)+(t(q2,npte)+(ht(q1)/dd(q2)))*w(q2,1:3);
        q2n=1+mod(q2,ne); %Next q2 index in cyclic order
        Z2 = m(q2n,1:3)+x(1)*X(q2n,1:3)+ ...%The first point on the next edge
            y(1)*Y(q2n,1:3)+(t(q2n,1)+(ht(q1)/dd(q2n)))*w(q2n,1:3);
        %Test code:
        Zdif = Z2-Z1;
        nZ =sqrt(Zdif*Zdif');
        Z0 = m(q2,1:3)+x(npte-1)*X(q2,1:3)+ ...%The last point on the current edge
            y(npte-1)*Y(q2,1:3)+(t(q2,npte-1)+(ht(q1)/dd(q2)))*w(q2,1:3);
        cnZ = sqrt((Z1-Z0)*(Z1-Z0)');
        %
        %
        Z(q2*npte+(q2-1)*nflp+1:q2*(npte+nflp),1:3) = ones(nflp,1)*Z1(1:3)+...
            ((1:nflp)'/(1+nflp))*Zdif(1:3);%Construct the straight segment joining the two curved ones
    end
    ntot = sz(Z); % This is the number of samples of the closed curve
%Before we Fourier transform we need to shift back to the plane <X,nv0>=0
 %   plot3(Z(1:ntot,1),Z(1:ntot,2),Z(1:ntot,3))
    Z(1:ntot,1:3) = Z(1:ntot,1:3)-ht(q1)*ones(ntot,1)*nv0(1:3);
    
    Zf=fft(Z); %This gives the Fourier transform of the columns of Z
    rho(q1) = (rs(q1)/ars(q1)); %This is the radius used in the harmonic extension
    nt2 = floor(ntot/2);
    nt3 = ntot-nt2;
    for q3 = 1:3 %We multiply the FT by the right power of rho to get a harmonic extension
        Zf(1:nt2,q3)=(rho(q1).^(0:nt2-1))'.*Zf(1:nt2,q3);
        Zf(nt2+1:ntot,q3)=(rho(q1).^(nt3:-1:1))'.*Zf(nt2+1:ntot,q3);
    end
    Zh = real(ifft(Zf));
    vrtxsamp{q1}(1:ntot,1:3)=Zh(1:ntot,1:3)+ht(q1)*ones(ntot,1)*nv0(1:3);
end
   
end
