function [Vout,Eout,Fout,lmn] = polyframe(V,E,F)

% This function inputs a geometric/combinatorial description of a
% polyhedron as: 
%
%     V = [[v1];[v2];...;[vnv]] is a matrix whose rows are
%           the vertices.
%     E = {[i11,i12],...,[ie1,ie2]} is a list of indices of
%           the endpoints of edges 
%     F = {[j1i,...,jin1],...,[jf1,...,jffnf]} is
%           list of the vertices that make up each face.
%
% It outputs data structures that facilitate that describe the local
% structure of the polyhedron in a neighborhood of each vertex, edge
% and face:
%
%     Vout ={{[i11,...,i1m],[nv1],{j11,...,j1m},...,{..}} For each
%          vertex we output the indices of the edges joined to it,
%          sorted in the correct cyclic order, an outward pointing
%          unit support vector, and the faces that define the edges in
%          the corresponding cyclic order:
%
%          e_{ikj} = f_{ikj} \cap f{ik(j+1)}.
%
%     Eout = {{[i11,i12],[mdpnt],[ne],[nv1,a1],[nv2,a2]},...{...}} For
%          each edge we output the indices of the faces that define
%          it, the midpoint of the edge, a normalized vector in the
%          direction of the edge, the equations of the faces defining
%          the faces:
%
%          <X,nvi> = ai.
%
%     Fout = {{[nv1],a},...,{[nvf],af}} nvj is the outer normal vector
%          to a face, and aj is the affine parameter:
%
%          The equation of the face is <X,nvj> = aj.
%


%
% set the intersection (numerical) tolerance
%
tol = 10^(-10);

%
% Find out how many of each object we have:
%
nv = sz(V);
ne = sz(E);
nf = sz(F);

%
% For later applications we find the centroid of the polyhedron
%
cv = sum(V)/nv;


%
% Process the faces
%
for j=1:nf

    nvf = sz(F{j});
    i1 = F{j}(1);
    i2 = F{j}(2);
    i3 = F{j}(3);
    v1 = V(i2,1:3)-V(i1,1:3);
    v2 = V(i2,1:3)-V(i3,1:3);
    
    %Compute the cross product to get the normal vector
    nvec = cp(v1,v2); 
    nrm = sqrt(nvec*nvec');
    nvec = nvec/nrm;
    
    % Make sure the orientation is correct by comparing <nv,cv> to
    % <nv,V(i1)>
    if cv*nvec'>V(i1,1:3)*nvec';
        nvec=-nvec;
    end
    ap = nvec*V(i1,1:3)';

    % Check that the vertices lie on the face and find the centroid.
    cent = [0,0,0];
    for j2 = 1:nvf
        if abs(V(F{j}(j2),1:3)*nvec'-ap)>tol
            error('The vertices in face %d do not lie on the face',j)
        end
        cent = cent+V(F{j}(j2),1:3);
    end

    %Outer normal vector
    Fout{j,1}=nvec; 

    %Affine parameter for face
    Fout{j,2}=ap; 

    %The centroid of the face
    Fout{j,3}(1:3)=cent(1:3)/nvf; 
end



%
%Process the edges
%
for j = 1:ne

    i1 = E{j}(1);
    i2 = E{j}(2);

    %The length of the edge
    len(j) = sqrt((V(i1,1:3)-V(i2,1:3))*(V(i1,1:3)-V(i2,1:3))'); 
    jf = 0;

    clear faces
    for j1 = 1:nf
        if asubb(i1,F{j1})+asubb(i2,F{j1}) == 2
            jf = jf +1;
            faces(jf) = j1;
        end
    end

    %Check that the edge lies on two faces
    if jf ~= 2
        error('Edge %d does not lie on two faces',j)
    end

    %Output the indices of the faces defining the edge.
    Eout{j,1} = faces; 

    %Output the midpoint of the edge.
    Eout{j,2} = (V(i1,1:3)+V(i2,1:3))/2; 
    ve = V(i1,1:3)-V(i2,1:3);

    %Output the normalized direction of the edge.
    Eout{j,3} = ve/sqrt(ve*ve'); 

    %Output the equations of the faces as [nvec, ap].
    Eout{j,4}= [Fout{faces(1),1},Fout{faces(1),2}]; 
    Eout{j,5}= [Fout{faces(2),1},Fout{faces(2),2}]; 

    %Output the indices of the vertices defining the edge.
    Eout{j,6} = [i1, i2]; 

    %This is the length of the shortest edge, which we use to scale all
    %subsequent lengths.
    lmn = min(len); 

end



%
%Process the vertices
%
for j = 1:nv
    clear edges ar ars ix Y
    % Find the other vertices joined to our vertex

    % At the end of the next loop this is the number of vertices
    % joined to the current vertex
    j1 = 0; 
    
    for j2 = 1:ne
        if E{j2}(1) == j
            j1 = j1+1;
            % Index of the other vertex
            verts(j1) = E{j2}(2);
            % Also keep track of the indices of the edges
            edges(j1) = j2; 
        elseif E{j2}(2) == j
            j1 = j1 + 1;
            % Index of the other vertex
            verts(j1) = E{j2}(1);
            % Also keep track of the indices of the edges
            edges(j1) = j2; 
        end
    end
    % Find an outer support vector at the vertex
    nv0=[0,0,0];

    for j3=1:j1
        % The indices of the next faces
        findx1 = Eout{edges(j3),1}(1); 
        % The indices of the next faces
        findx2 = Eout{edges(j3),1}(2); 
        %Add the faces' support vectors
        nv0 = nv0 + Fout{findx1,1}(1:3) + Fout{findx2,1}(1:3); 
    end

    nv0 = nv0/sqrt(nv0*nv0');

    % Output the support vector at the vertex.
    Vout{j,2} = nv0; 
    
    % Now we need to find the correct cyclic order for these vertices.
    for j3 = 1:j1
        vn = (V(verts(j3),1:3)-V(j,1:3));
        vn = vn/sqrt(vn*vn');
        Y(j3,1:3) = vn - (vn*nv0')*nv0;
        Y(j3,1:3) = Y(j3,1:3)/sqrt(Y(j3,1:3)*Y(j3,1:3)');
    end
    
    % Choose and o/n basis for the plane <X,nv0>=0
    e1 = Y(1,1:3);
    e2 = cp(nv0,e1);
    for j3 = 1:j1
        ar(j3) = angle(e1*Y(j3,1:3)'+i*(e2*Y(j3,1:3)'));
    end
    [ars,ix] = sort(ar);

    % Output the sorted indices of the vertices joined to our vertex.
    Vout{j,1}=verts(ix); 
    
    % Now sort the indices of faces to be consistent with sorting of the
    % edges.
    for j3 = 1:j1
        i1 = Eout{edges(ix(j3)),1}(1);
        i2 = Eout{edges(ix(j3)),1}(2);
        j3n = 1+mod(j3,j1);
        i3 = Eout{edges(ix(j3n)),1}(1);
        i4 = Eout{edges(ix(j3n)),1}(2);
        if i1 == i3 | i1 == i4
            Vout{j,3}(j3) = i2;
        elseif i2 == i3 | i2 == i4
            Vout{j,3}(j3) = i1;
        else error('Edge %d on vertex %d is out of order', j3,j)
        end
    end
    
end
 




