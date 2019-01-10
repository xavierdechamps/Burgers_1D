clear -all
clc
close all
format long
warning off

Length_domain    = 2*pi  % 1 or 2 for the sinus wave
Viscosity        = 0.0075 % 0.001 or higher for the sinus wave
Number_elements  = 128
Subgrid_constant = 0.;
Time_total       = 40 % 210 % 1 for the sinus wave
Time_steps       = 4000
Ninterpolation   = 50;
Name_output      = 'test' %FD_spectral_like_optimal_7.5e-3_2048
Name_spectrum_ref= 'Results/spectrum_2048_35e-4_new.txt';

addpath('./src/');

%FE_LagrangeP1(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

%FE_LagrangeP3(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);

%FE_HermiteP3 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);

%FE_HermiteP5 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Ninterpolation,Name_output,Name_spectrum_ref);

%FD_conservative_order2(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

%FD_conservative_order4(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

FD_dissipative_order2(Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref);

%FD_upwind_order1 (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref) ;

%FD_compact_spectral (Number_elements,Viscosity,Subgrid_constant,Length_domain,Time_total,Time_steps,Name_output,Name_spectrum_ref) ;
