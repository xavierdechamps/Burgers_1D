Solve the 1D random forced viscous Burgers equation with high order finite element and finite difference methods.

The main driver of the code is the file Main.m

The implemented finite element methods are:
  * linear Lagrange element
  * cubic Lagrange element
  * cubic Hermite element
  * 5th order Hermite element
  
The implemented finite difference schemes are
  * energy dissipative second order centered scheme
  * energy conservative second order centered scheme
  * energy conservative fourth order centered scheme
  * energy conservative compact schemes with spectral-like resolution
  * non-linear discretization of the convective term (slope-limiters)
  
The implemented finite difference slope limiters are chosen from the article "On the spectral and conservation properties of nonlinear discretization operators" by D. Fauconnier and E. Dick. These are the
  * central DFD 
  * 1st, 2nd and 3rd order upwind discretization (UP1, UP2 and UP3)
  * 5th order Weighted Essentially Non-Oscillatory (WENO)
         
The energy spectrum is computed and compared with a reference spectrum from a pseudo-spectral code.