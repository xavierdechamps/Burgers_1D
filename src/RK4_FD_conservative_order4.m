function y = RK4_FD_conservative_order4 (u,deltat,N,K,F,h,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with an energy conservative Hs scheme of 
% order 4 for the convective term
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

  ind = zeros(N,5);
  ind(:,1) = -1:N-2; 
  ind(:,2) = 0:N-1; 
  ind(:,3) = 1:N;
  ind(:,4) = 2:N+1;
  ind(:,5) = 3:N+2;
  ind(1,:)   = [N-1 N 1 2 3]; 
  ind(2,:)   = [N   1 2 3 4]; 
  ind(N,:)   = [N-2 N-1 N   1 2]; 
  ind(N-1,:) = [N-3 N-2 N-1 N 1]; 

  sixth = 1/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C = get_nonlinear_term(Un,ind,h,N);
  k1 = ( K - C ) * Un ;
  Un2 = Un + deltat*0.5*k1 ;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C = get_nonlinear_term(Un2,ind,h,N);
  k2 = ( K - C ) * Un2;
  Un3 = Un + deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C = get_nonlinear_term(Un3,ind,h,N);
  k3 = ( K - C ) * Un3;
  Un4 = Un + deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C = get_nonlinear_term(Un4,ind,h,N);
  k4 = ( K - C ) * Un4 ; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y = Un + deltat*sixth*(k1 + 2*k2 + 2*k3 + k4 ) + F ;
  
end

function vecC = get_nonlinear_term(Un,ind,h,N)
  thirtysixth_over_h = 1/(36*h);
   
% advective - on diagonal
  diag0  =  Un(ind(:,1)) - 8*Un(ind(:,2)) + 8*Un(ind(:,4)) - Un(ind(:,5)) ;
% dissipative - off diagonal
  diagm2 =      Un(ind(:,1));
  diagm1 = -8 * Un(ind(:,2));
  diagp1 =  8 * Un(ind(:,4));
  diagp2 =    - Un(ind(:,5));

  C0     = diag(sparse(diag0)); Cp1 = diag(sparse(diagp1),1); Cm1 = diag(sparse(diagm1),-1); Cp2 = diag(sparse(diagp2),2); Cm2 = diag(sparse(diagm2),-2);
  C      = C0 + Cp1(1:N,1:N) + Cm1(2:N+1,2:N+1) + Cp2(1:N,1:N) + Cm2(3:N+2,3:N+2);
  C(N,1) = diagp1(N);   C(N,2)   = diagp2(N);    C(N-1,1) = diagp2(N-1);
  C(1,N) = diagm1(1);   C(1,N-1) = diagm2(1);    C(2,N)   = diagm2(2);
  
  vecC = C * thirtysixth_over_h;
endfunction
