function y = RK4_FD_DFD (u,deltat,N,K,F,h,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with a nonlinear discretization DFD for the convective term
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

  ind = zeros(N,7);
  ind(:,1) = -2:N-3; 
  ind(:,2) = -1:N-2; 
  ind(:,3) = 0:N-1;
  ind(:,4) = 1:N;
  ind(:,5) = 2:N+1;
  ind(:,6) = 3:N+2;
  ind(:,7) = 4:N+3;
  ind(1,:) = [N-2 N-1 N 1 2 3 4]; 
  ind(2,:) = [N-1 N   1 2 3 4 5]; 
  ind(3,:) = [N   1   2 3 4 5 6]; 
  ind(N-2,:) = [N-5 N-4 N-3 N-2 N-1 N 1];
  ind(N-1,:) = [N-4 N-3 N-2 N-1 N   1 2];
  ind(N,:)   = [N-3 N-2 N-1 N   1   2 3];

  f = 0.2;
  one_over_h = 0.25/(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C = - get_nonlinear_term(Un,ind,f) * one_over_h ;
  k1  = K * Un + C;
  Un2 = Un + deltat * 0.5 * k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C = - get_nonlinear_term(Un2,ind,f) * one_over_h ;
  k2  = K * Un2 + C;
  Un3 = Un + deltat * 0.5 * k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C = - get_nonlinear_term(Un3,ind,f) * one_over_h ;
  k3  = K * Un3 + C;
  Un4 = Un + deltat * k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C = - get_nonlinear_term(Un4,ind,f) * one_over_h ;
  k4 = K * Un4 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y  = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 +F;

end

function vecC = get_nonlinear_term(Un,ind,f)
  cip1 = - 1 / ( 6 * (1 + f * max( min( ( Un(ind(:,7))-4*Un(ind(:,6))+6*Un(ind(:,5))-4*Un(ind(:,4))+Un(ind(:,3)) )./( Un(ind(:,6))-2*Un(ind(:,5))+Un(ind(:,4)) ) , 0 ) , -3 ) ) );
  cim1 = - 1 / ( 6 * (1 + f * max( min( ( Un(ind(:,5))-4*Un(ind(:,4))+6*Un(ind(:,3))-4*Un(ind(:,2))+Un(ind(:,1)) )./( Un(ind(:,4))-2*Un(ind(:,3))+Un(ind(:,2)) ) , 0 ) , -3 ) ) );
    
  vecC = Un(ind(:,5)).^2 + ...
         cip1'.*( Un(ind(:,6)).^2 -2*Un(ind(:,5)).^2 + Un(ind(:,4)).^2 ) - ...
         Un(ind(:,3)).^2 - ...
         cim1'.*( Un(ind(:,4)).^2 -2*Un(ind(:,3)).^2 + Un(ind(:,2)).^2 )  ;

endfunction
