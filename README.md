Solve the 1D random forced viscous Burgers equation with high order finite element and finite difference methods.
Direct Numerical Simulation and Large-Eddy Simulation are possible. The turbulence model implemented for LES is the eddy viscosity Smagorinsky model.
Both a constant Smagorinsky model and a dynamic Smagorinsky model are implemented.

The main driver of the code is the file Main.m

The implemented finite element methods are:
  * continuous and discontinuous linear Lagrange element
  * continuous 3rd order Lagrange element
  * continuous and discontinuous 3rd order Hermite element
  * continuous 5th order Hermite element
  
The implemented finite difference schemes are
  * energy dissipative 2nd order centered scheme
  * energy conservative 2nd order centered scheme
  * energy conservative 4th order centered scheme
  * energy conservative compact schemes with spectral-like resolution
  * non-linear discretization of the convective term (slope-limiters)
  
The implemented finite difference slope-limiters are chosen from the article "On the spectral and conservation properties 
of nonlinear discretization operators" by D. Fauconnier and E. Dick. These are the
  * central Dynamic Finite Difference (DFD) 
  * 1st, 2nd and 3rd order upwind discretization (UP1, UP2 and UP3)
  * Total Variation Diminishing (TVD) scheme
  
The non-linear discretization is based on the skew-symmetric form of the convective term. 
         
The energy spectrum is computed and compared with a reference spectrum from a pseudo-spectral code.