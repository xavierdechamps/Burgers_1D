function [y,energy] = RK4_FD_compact_spectral (u,deltat,N,mat_deriv1,mat_deriv2,F,h,constant_sub,factors_deriv)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with compact finite difference schemes
% 
% The equation to be solved is 
%                  du
%                  -- = f(u,t)
%                  dt
% The explicit 4 steps Runge-Kutta scheme is 
% U(n+1) = U(n) + 1/6(k1 + 2k2 + 2k3 +k4)
% where k1 = f(U(n),t)
%       k2 = f(U(n) + (deltat/2)k1,t + deltat/2)
%       k3 = f(U(n) + (deltat/2)k2,t + deltat/2)
%       k4 = f(U(n) + deltat.k3,t + deltat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
  Un = u;
    
  ind = zeros(N,7);
  ind(:,1) = -2:N-3; 
  ind(:,2) = -1:N-2; 
  ind(:,3) = 0:N-1; 
  ind(:,4) = 1:N;
  ind(:,5) = 2:N+1;
  ind(:,6) = 3:N+2;
  ind(:,7) = 4:N+3;
  ind(N-2,:) = [N-5 N-4 N-3 N-2 N-1 N 1]; 
  ind(N-1,:) = [N-4 N-3 N-2 N-1 N   1 2]; 
  ind(N,:)   = [N-3 N-2 N-1 N   1   2 3]; 
  ind(1,:)   = [N-2 N-1 N   1   2   3 4]; 
  ind(2,:)   = [N-1 N   1   2   3   4 5]; 
  ind(3,:)   = [N   1   2   3   4   5 6]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% Skew-symmetric form of the convective term
% First consider the divergence form
  vec_deriv1 = get_first_derivative_div(Un,ind,factors_deriv);
  nonlin_div = mat_deriv1\vec_deriv1;

% Then the advective form
  deriv1_u = get_first_derivative_adv(Un,ind,factors_deriv);
  nonlin_adv = Un.* ( mat_deriv1\deriv1_u );
  
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un,ind,factors_deriv);
  
  k1 = -(nonlin_adv + 2*nonlin_div)/3 + mat_deriv2\vec_deriv2;
  Un2 = Un+deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
% Skew-symmetric form of the convective term
% First consider the divergence form
  vec_deriv1 = get_first_derivative_div(Un2,ind,factors_deriv);
  nonlin_div = mat_deriv1\vec_deriv1;
  
% Then the advective form
  deriv1_u = get_first_derivative_adv(Un2,ind,factors_deriv);
  nonlin_adv = Un2.* ( mat_deriv1\deriv1_u );
    
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un2,ind,factors_deriv);
  
  k2 = -(nonlin_adv + 2*nonlin_div)/3 + mat_deriv2\vec_deriv2;
  Un3 = Un+deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
% Skew-symmetric form of the convective term
% First consider the divergence form
  vec_deriv1 = get_first_derivative_div(Un3,ind,factors_deriv);
  nonlin_div = mat_deriv1\vec_deriv1;  
  
% Then the advective form
  deriv1_u = get_first_derivative_adv(Un3,ind,factors_deriv);
  nonlin_adv = Un3.* ( mat_deriv1\deriv1_u );
    
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un3,ind,factors_deriv);

  k3 = -(nonlin_adv + 2*nonlin_div)/3 + mat_deriv2\vec_deriv2;
  Un4 = Un+deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
% Skew-symmetric form of the convective term
% First consider the divergence form
  vec_deriv1 = get_first_derivative_div(Un4,ind,factors_deriv);
  nonlin_div = mat_deriv1\vec_deriv1;
  
% Then the advective form
  deriv1_u = get_first_derivative_adv(Un4,ind,factors_deriv);
  nonlin_adv = Un4.* ( mat_deriv1\deriv1_u );
  
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un4,ind,factors_deriv);
    
  k4 = -(nonlin_adv + 2*nonlin_div)/3 + mat_deriv2\vec_deriv2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  y = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 + F;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%    energy = 0;
    energy = get_energy(y,ind,mat_deriv1,mat_deriv2,factors_deriv,h,F);

endfunction


function dfdx = get_first_derivative_div(Un,ind,factors_deriv)
% Compute the first spatial derivative with a compact finite difference scheme, in the divergence form
  dfdx = 0.5 * factors_deriv(5,1) * ( Un(ind(:,7)).*Un(ind(:,7)) - Un(ind(:,1)).*Un(ind(:,1)) ) + ...
         0.5 * factors_deriv(4,1) * ( Un(ind(:,6)).*Un(ind(:,6)) - Un(ind(:,2)).*Un(ind(:,2)) ) + ...
         0.5 * factors_deriv(3,1) * ( Un(ind(:,5)).*Un(ind(:,5)) - Un(ind(:,3)).*Un(ind(:,3)) ) ;
endfunction


function dfdx = get_first_derivative_adv(Un,ind,factors_deriv)
% Compute the first spatial derivative with a compact finite difference scheme, in the advective form
  dfdx = factors_deriv(5,1) * ( Un(ind(:,7)) - Un(ind(:,1)) ) + ...
         factors_deriv(4,1) * ( Un(ind(:,6)) - Un(ind(:,2)) ) + ...
         factors_deriv(3,1) * ( Un(ind(:,5)) - Un(ind(:,3)) ) ;
endfunction


function d2fdx2 = get_second_derivative(Un,ind,factors_deriv)
% Compute the second spatial derivative with a compact finite difference scheme
  d2fdx2 = factors_deriv(5,2) * ( Un(ind(:,7)) -2*Un(ind(:,4)) + Un(ind(:,1)) ) + ...
           factors_deriv(4,2) * ( Un(ind(:,6)) -2*Un(ind(:,4)) + Un(ind(:,2)) ) + ...
           factors_deriv(3,2) * ( Un(ind(:,5)) -2*Un(ind(:,4)) + Un(ind(:,3)) ) ;  
endfunction

function energy = get_energy(Un,ind,mat_deriv1,mat_deriv2,factors_deriv,h,F)
% Compute the numerical energy produced by the spatial discretization of the 
% convective term in the skew-symmetric form
  energy = 0;
  
% Divergence form
  vec_deriv1 = get_first_derivative_div(Un,ind,factors_deriv);
  deriv1_u = mat_deriv1\vec_deriv1;
  energy1 = h * Un' * deriv1_u;
    
% Advective form
  vec_deriv1 = get_first_derivative_adv(Un,ind,factors_deriv);
  deriv1_u = mat_deriv1\vec_deriv1;
  energy2 = h * (Un .* Un)' * deriv1_u; 

% Sum between advective and divergent forms for nonlinear term
  energy = energy - (2*energy1 + energy2)/3 ;
    
% % Contribution for viscous dissipation
%  vec_deriv2 = get_second_derivative(Un,ind,factors_deriv) ;
%  deriv2_u = mat_deriv2\vec_deriv2;
%  energy = energy + h * Un' * deriv2_u;
     
% % Contribution of forcing term
%  energy = energy + h * Un' * F;
endfunction