function y = RK4_FD_conservative_order2 (u,deltat,N,K,F,h,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with an energy conservative Hs scheme of 
% order 2 for the convective term
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
	ind(1,:)   = [N 1   2]; 
	ind(N,:) = [N-1 N 1];

	sixth = 1/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C = get_nonlinear_term(Un,ind,constant_sub,h,N);
	k1=(K+C)*Un;
	Un2 = Un+deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C = get_nonlinear_term(Un2,ind,constant_sub,h,N);
	k2=(K+C)*Un2;
	Un3 = Un+deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C = get_nonlinear_term(Un3,ind,constant_sub,h,N);
	k3=(K+C)*Un3;
	Un4 = Un+deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C = get_nonlinear_term(Un4,ind,constant_sub,h,N);
	k4=(K+C)*Un4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y = Un + deltat*sixth*(k1 + 2*k2 + 2*k3 +k4 ) +F;
	
end

function vecC = get_nonlinear_term(Un,ind,constant_sub,h,N)
	sixth = 1/6;
	one_over_h = 1/h;
	one_over_6h = sixth*one_over_h;
  
% non-linear terms
	diag0  = one_over_6h * ( Un(ind(:,1),1) - Un(ind(:,3),1) );
	diagm1 = diag0 ;
  diagp1 = diag0 ;

% subgrid terms
	diagm1 = diagm1 + (constant_sub * one_over_h) *   abs(Un(ind(:,2),1) - Un(ind(:,1),1));
	diag0  = diag0  - (constant_sub * one_over_h) * ( abs(Un(ind(:,2),1) - Un(ind(:,1),1)) + abs(Un(ind(:,3),1) - Un(ind(:,2),1)) );
	diagp1 = diagp1 + (constant_sub * one_over_h) *   abs(Un(ind(:,3),1) - Un(ind(:,2),1));

	C0  = diag(sparse(diag0));
  Cp1 = diag(sparse(diagp1),1);
  Cm1 = diag(sparse(diagm1),-1);
  
	vecC      = C0 + Cp1(1:N,1:N) + Cm1(2:N+1,2:N+1);
	vecC(N,1) = diagp1(N,1);
	vecC(1,N) = diagm1(1,1);
endfunction
