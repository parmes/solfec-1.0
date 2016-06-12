# Solfec
Parallel Non-Smooth Contact Dynamics using MPI, C and Python.

Solfec implements the Non-Smooth Contact Dynamics Method [1,2] using MPI, C, Python and several 3rd party
codes written in C/C++/Fortran. It includes mesh, convex polyhedra, sphere and ellipsoid based shapes,
linear elastic first order finite elements, pseudo-rigid and rigid kinematics, velocity based Signorini-Coulomb
contact/impact law, and a parallel time stepping combined with a simple dynamic load balancing. Solfec
has been developed as a part of research [3], and it is currently developed and maintained with support
from the civil nuclear context int the UK.

[1] J. J. Moreau. Numerical aspects of the sweeping process. CMAME, 177(3-4):329–349, 1999.
[2] M. Jean. The non-smooth contact dynamics method. CMAME, 177(3-4):235–257, 1999.
[3] T. Koziara, N. Bićanić. A distributed memory parallel multibody Contact Dynamics code. IJNME, 87(1-5):437-456, 2011.
