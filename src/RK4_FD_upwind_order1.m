function y = RK4_FD_upwind_order1 (u,deltat,N,K,F,h,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with a nonlinear discretization first order UP1 scheme for the convective term
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
  Un=u;

  ind = zeros(N,3);
  ind(:,1) = 0:N-1; 
  ind(:,2) = 1:N; 
  ind(:,3) = 2:N+1;
  ind(1,:) = [N   1 2]; 
  ind(N,:) = [N-1 N 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C = - get_nonlinear_term(Un,ind,h) ;
  k1  = K * Un + C;
  Un2 = Un + deltat * 0.5 * k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C = - get_nonlinear_term(Un2,ind,h) ;
  k2  = K * Un2 + C;
  Un3 = Un + deltat * 0.5 * k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C = - get_nonlinear_term(Un3,ind,h) ;
  k3  = K * Un3 + C;
  Un4 = Un + deltat * k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C = - get_nonlinear_term(Un4,ind,h) ;
  k4 = K * Un4 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y  = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 +F;
  
end

function vecC = get_nonlinear_term(Un,ind,h)
% i = i-1  i  i+1
%      1   2   3  

  vpi = max( Un(ind(:,2),1) ,0); % v+ at node i
  vmi = min( Un(ind(:,2),1) ,0); % v- at node i
 
  vpip1 = max( ( Un(ind(:,2),1) + Un(ind(:,3),1) ) *0.5   ,0); % v+ at node i+1/2
  vpim1 = max( ( Un(ind(:,2),1) + Un(ind(:,1),1) ) *0.5   ,0); % v+ at node i-1/2
  vmip1 = min( ( Un(ind(:,2),1) + Un(ind(:,3),1) ) *0.5   ,0); % v- at node i+1/2
  vmim1 = min( ( Un(ind(:,2),1) + Un(ind(:,1),1) ) *0.5   ,0); % v- at node i-1/2
  
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = ( ( vpi + vpip1 - vmi - vmim1 ) .* Un(ind(:,2),1) + ...
           ( vmi + vmip1 ) .* Un(ind(:,3),1) - ...
           ( vpi + vpim1 ) .* Un(ind(:,1),1) ) / ( 3. * h ) ;
 
% Contains only the divergence form -> not energy conservative -> not good!
%  vecC = ( (vpip1-vmim1).*Un(ind(:,2),1) + ...
%           vmip1.*Un(ind(:,3),1) - ...
%           vpim1.*Un(ind(:,1),1) ) * 0.5 / h ;
endfunction
