clear -all
clc
close all
format long
warning off

Length_domain     = 2*pi   % Total length of the spatial domain (put 1 or 2 for the sinus wave as initial condition)
Viscosity         = 0.0035 % Kinematic viscosity
Number_elements   = 128     % Number of elements (for FE) / nodes (for FD)
Subgrid_constant  = 0.5   % Subgrid terms are implemented for some discretizations
Filter            = 1      % Dynamic Smagorinsky model: type of filter (0 = deactivated -> use the constant Subgrid_constant)
                                   % 1 = 3-nodes binomial filter
                                   % 2 = 5-nodes binomial filter
                                   % 3 = 7-nodes binomial filter
                                   % 4 = 9-nodes binomial filter
Order_Viscous     = 6       % = 2, 4, 6 : Order of the discretization for second derivative (used in Hc4 and in nonlinear schemes)
Time_total        = 410     % Time of the simulation (1 for the sinus wave)
Time_steps        = 40000   % Number of time steps, increment in time is thus equal to Time_total/Time_steps
Ninterpolation    = 20;    % The high-order finite elements show oscillations within the elements,
                           % This parameter indicates how many additionnal points are interpolated for each element
Name_output       = 'test' % Name of the output file
Name_spectrum_ref = 'Results/spectrum_2048_35e-4_new.txt'; % Name of the reference turbulence spectrum calculated by a pseudo-spectral method
% Results/spectrum_2048_35e-4_new.txt has 1024 modes and has been calculated for a kinematic viscosity = 0.0035
% Results/spectrum_2048_75e-4_new.txt has 1024 modes and has been calculated for a kinematic viscosity = 0.0075

% Method = 1 : finite element (FE) linear Lagrange
% Method = 2 : finite element (FE) cubic Lagrange
% Method = 3 : finite element (FE) cubic Hermite
% Method = 4 : finite element (FE) 5th order Hermite
% Method = 5 : finite difference (FD) energy conservative order 2 for convective term
% Method = 6 : finite difference (FD) energy conservative order 4 for convective term
% Method = 7 : finite difference (FD) energy dissipative order 2 for convective term
% Method = 8 : finite difference (FD) compact spectral-like resolution
% Method = 9 : finite difference (FD) non-linear discretization of the convective term
method = 1 ;
submethod = 1 ; % Used by the compact FD, the nonlinear FD schemes and cubic Lagrange + Hermite FE

addpath('./src/');

switch method
  case 1
     FE_LagrangeP1(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);
    
  case 2
     FE_LagrangeP3(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref,submethod);
  
  case 3
     FE_HermiteH3 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref,submethod);

  case 4
     FE_HermiteH5 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref,submethod);

  case 5
     FD_conservative_order2(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

  case 6
     FD_conservative_order4(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref,Order_Viscous);

  case 7
     FD_dissipative_order2(Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

  case 8
     FD_compact_spectral (Number_elements,Viscosity,Subgrid_constant,Filter,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref,submethod,Order_Viscous) ;

  case 9
     FD_nonlinear_schemes(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref,submethod,Order_Viscous) ;

  otherwise
     disp('Method is not implemented yet')
     
end