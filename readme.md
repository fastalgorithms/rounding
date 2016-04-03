## Polyhedral smoothing in MATLAB

![alt text](pyramid.png "Rounded pyramid")

The code is divided into three main applications:

### Smoothing 1-dimensional objects in the plane

`smthabv(h,k,n)` This program returns $2n$ equispaced samples of the convolution $|x|*(1-(x/h)^2)^k$ on $[-1,1]$

`smthpoly(v,h,k,m)` This program smooths a polygon with vertices listed in the variable `v`. The other parameters `h,k,m` specify the smoothness of the resulting smoothed vertices and the number of sample points therein.

These programs are called by `polygon_tst.m` which uses these routines to smooth regular polygons and plot the results.

### Smoothing 2-dimensional objects in 3-space.

a. polyframe(V,E,F) takes a polyhedron with vertices in V and edges
   and faces described in E and F respectively and returns data
   structures describing the local geometry of vertices, edges and
   faces. For example, for each vertex we get a list of the other
   vertices joined to the given vertex by an edge, and a choice of
   outer point unit vector at the vertex, and the faces that define
   each edge. This facilitates the computations used to smooth the
   polyhedron.
      b. smthvrtx(i,V,E,Vout,Fout,he,ke,hv,kv,me,mv) This program
   smooths the polyhedron near to the vertex i. As before (V,E,Vout,
   Fout) are geometric and combinatorial information about the
   polyhedron and (he,ke,hv,kv,me,mv) are parameters that specify the
   smoothing of the edges, (he,ke), and vertices, (hv,kv), as well as
   the number of samples returned on the smoothed objects (me,mv).
      c. smthedge(v0,v1,nv0,nv1,mv0,mv1,lmin0,lmin1,h,k,ne,nh) is a
   program that smooths the edges of polyhedra. This program assumes
   that the vertices have already been smoothed, and the variables
   (lmin0,lmin1) describe the size of the region along edge modified
   in the vertex smoothing step.

These programs are called by smth_polyhedron_tst.m, which contains the
data needed to define the platonic solids (but not the dodecahedron).

### Plotting

a. pltsmvx plots samples of a smoothed vertex
b. pltsmed plots samples of a smoothed edge
c. pltpoly plots the edges and vertices of a polyhedron

These programs are used by smth_polyhedron_tst.m to display the output.
