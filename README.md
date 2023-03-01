# Bachelor-Thesis-Auxiliary-Code
I decided to upload all the code I've been writing for my bachelor thesis to have it accessible from anywhere and stored somewhere safe.

This probably won't be of much use for anyone besides me and next generations of students programming numerics for [Template Numerical Library](https://gitlab.com/tnl-project/tnl).

The code is a mixture of my own contributions, samples from TNL tutorials and the recommended *Makefile* and *config.mk* provided at TNL's GitLab

Short description of files follows:

**aux.cpp**: me getting to know how to initialize containers. The templates can get tricky though, so it's useful to have a backup for later work.

**iter.cpp**: iterates over all the cells and faces of a mesh and returns their respective n-dimensional Lebesgue measures, as well as sums of those, separately for boundary and interior entities.

**normals.cpp**: iterating over all the cells and then over each cell's faces, this code returns *unit outward normal vectors* for each cell's face in the form of *TNL::Containers::Vector* of smaller *TNL::Containers::Vector* entities, which contain the normals of type *PointType*, which is again a *TNL::Containers::Vector*. The word vector is a little overloaded here.

**numgrad.cpp**: the backbone of it all, this is the numerical scheme which computes gradient inside a cell. Based on Generalized Stokes' Theorem, it translates a volume integral of *grad(f)* to a surface integral of just *f*. This approximation converges to the analytical (real) gradient of *f*. The convergence is tested on analytical functions, whose analytical gradient is calculated manually in the function *angrad*.

**autodiff.ipynb**: tests the handcalculated derivatives in comparison with the results of differentiating with the python library mygrad
