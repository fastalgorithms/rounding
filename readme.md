## Polyhedral smoothing in MATLAB

![alt text](icos.png "Rounded icosahedron vertex")

The code is divided into three main applications:

### Smoothing 1-dimensional objects in the plane

The driver file **polygon_test.m** constructs rounding of regular
n-gons. It is easily adaptable to other polygons, only vertex
locations are needed. The resulting polygons are automatically
plotted.

The main routines used are `smthabv(h,k,n)` and `smthpoly(v,h,k,m)`.

The file [smthabv.m](src/smthabv.m) contains a function which
returns 2n equispaced samples of the convolution of a compactly
supported kernel ( 1-x^2 )^k on [-1,1] with the function abs(x).

The file [smthpoly.m](src/smthpoly.m) smooths a polygon with vertices
listed in the variable `v`. The other parameters `h,k,m` specify the
smoothness of the resulting smoothed vertices and the number of sample
points therein.


### Smoothing 2-dimensional objects in 3-space.

The driver file **polyhedron_test.m** constructs rounding of the
Platonic solids (except the dodecahedron).
It is easily adaptable to other polyhedrons, only vertex, edges, and
face information is needed. A carcass of the resulting polyhedron is
plotted automatically.

The main files used are [polyframe.m](src/polyframe.m),
[smthvrtx.m](src/smthvrtx.m),
and [smthedge.m](src/smthedge.m).

- `polyframe(V,E,F)` takes a polyhedron with vertices in `V`, edges
   and faces described in `E` and `F`, respectively, and returns data
   structures describing the local geometry of vertices, edges, and
   faces. For example, for each vertex a list of the other vertices
   joined to the given vertex by an edge is computed, along with an
   outward unit vector at the vertex, and faces that define each
   edge. This facilitates the computations used to smooth the
   polyhedron.
- `smthvrtx(i,V,E,Vout,Fout,he,ke,hv,kv,me,mv)` smooths the polyhedron
   near the vertex `i`. As before `V,E,Vout,Fout` are geometric and
   combinatorial information about the polyhedron and
   `he,ke,hv,kv,me,mv` are parameters that specify the smoothing of
   the edges `he,ke`, vertices `hv,kv`. The number of samples returned
   on the smoothed objects are set by `me,mv`.
- `smthedge(v0,v1,nv0,nv1,mv0,mv1,lmin0,lmin1,h,k,ne,nh)` is a
   subroutine that smooths the edges of polyhedra. This program
   assumes that the vertices have already been smoothed, and that the
   variables `lmin0,lmin1` describe the size of the region along the
   edge modified in the vertex smoothing step.

### Plotting

There are several utility plotting rountines:
- [pltsmvx.m](src/pltsmvx.m) plots samples of a smoothed vertex
- [pltsmed.m](src/pltsmed.m) plots samples of a smoothed edge
- [pltpoly.m](src/pltpoly.m) plots the edges and vertices of a polyhedron

