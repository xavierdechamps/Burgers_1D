clear -all
clc
close all
format long
warning off

Length_domain     = 2*pi   % Total length of the spatial domain (1 or 2 for the sinus wave)
Viscosity         = 0.0075 % Kinematic viscosity
Number_elements   = 64     % Number of elements (for FE) / nodes (for FD)
Subgrid_constant  = 0.2   % Subgrid terms are implemented for some discretizations
Filter            = 1      % Dynamic Smagorinsky model: type of filter (0 = deactivated -> use the constant Subgrid_constant)
                                   % 1 = 3 nodes binomial filter
                                   % 2 = 5 nodes binomial filter
                                   % 3 = 7 nodes binomial filter
                                   % 4 = 9 nodes binomial filter
Time_total        = 410     % Time of the simulation (1 for the sinus wave)
Time_steps        = 20000   % Number of time steps, increment in time is thus equal to Time_total/Time_steps
Ninterpolation    = 50;    % The high-order finite elements show oscillations within the elements, this parameter indicates how many additionnal points are interpolated for each element
Name_output       = 'essai' % Name of the output file
Name_spectrum_ref = 'Results/spectrum_2048_75e-4_new.txt'; % Name of the reference turbulence spectrum calculated by a pseudo-spectral method
% Results/spectrum_2048_35e-4_new.txt has 1024 modes and has been calculated for a kinematic viscosity = 0.0035
% Results/spectrum_2048_75e-4_new.txt has 1024 modes and has been calculated for a kinematic viscosity = 0.0075

% Method = 1 : finite element linear Lagrange
% Method = 2 : finite element cubic Lagrange
% Method = 3 : finite element cubic Hermite
% Method = 4 : finite element 5th order Hermite
% Method = 5 : finite difference energy conservative order 2 for convective term
% Method = 6 : finite difference energy conservative order 4 for convective term
% Method = 7 : finite difference energy dissipative order 2 for convective term
% Method = 8 : finite difference compact spectral-like resolution
% Method = 9 : finite difference non-linear discretization of the convective term
method = 7 ;
submethod = 1 ; % Used by the compact finite difference and by the nonlinear finite difference schemes

addpath('./src/');

switch method
  case 1
     FE_LagrangeP1(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);
    
  case 2
     FE_LagrangeP3(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);
  
  case 3
     FE_HermiteP3 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);

  case 4
     FE_HermiteP5 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);

  case 5
     FD_conservative_order2(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

  case 6
     FD_conservative_order4(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

  case 7
     FD_dissipative_order2(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

  case 8
     FD_compact_spectral (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref,submethod) ;

  case 9
     FD_nonlinear_schemes(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref) ;

  otherwise
     disp('Method is not implemented yet')
     
end