function y = RK4_FE_Lagrangep1 (u,deltat,N,M,k,nu,F,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with linear Lagrange elements
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

% Initialization of the non-linear term
  Cj=zeros(N,1);
  sixth = 1/6;
  half = 0.5;

  ind = zeros(N,3);
  ind(:,1) = 0:N-1; 
  ind(:,2) = 1:N; 
  ind(:,3) = 2:N+1;
  ind(1,:) = [N   1 2]; 
  ind(N,:) = [N-1 N 1];

  Un = u;
% nonlinear terms
  Cj = sixth * ( Un(ind(:,3)) .* ( Un(ind(:,2)) + Un(ind(:,3)) ) - Un(ind(:,1)) .* ( Un(ind(:,2)) + Un(ind(:,1)) ) );
% subgrid terms
  Cj = Cj - constant_sub * ( abs(Un(ind(:,3))-Un(ind(:,2))).*(Un(ind(:,3))-Un(ind(:,2))) - abs(Un(ind(:,2))-Un(ind(:,1))).*(Un(ind(:,2))-Un(ind(:,1))) );
  k1 = - M \ (nu*k*Un + Cj);

  Un2 = Un + deltat*half*k1;
% nonlinear terms
  Cj  = sixth * ( Un2(ind(:,3)) .* ( Un2(ind(:,2)) + Un2(ind(:,3)) ) - Un2(ind(:,1)) .* ( Un2(ind(:,2)) + Un2(ind(:,1)) ) );
% subgrid terms
  Cj = Cj - constant_sub * ( abs(Un2(ind(:,3))-Un2(ind(:,2))).*(Un2(ind(:,3))-Un2(ind(:,2))) - abs(Un2(ind(:,2))-Un2(ind(:,1))).*(Un2(ind(:,2))-Un2(ind(:,1))) );
  k2 = - M \ (nu*k*Un2 + Cj);

  Un3 = Un + deltat*half*k2;
% nonlinear terms
  Cj  = sixth * ( Un3(ind(:,3)) .* ( Un3(ind(:,2)) + Un3(ind(:,3)) ) - Un3(ind(:,1)) .* ( Un3(ind(:,2)) + Un3(ind(:,1)) ) );
% subgrid terms
  Cj = Cj - constant_sub * ( abs(Un3(ind(:,3))-Un3(ind(:,2))).*(Un3(ind(:,3))-Un3(ind(:,2))) - abs(Un3(ind(:,2))-Un3(ind(:,1))).*(Un3(ind(:,2))-Un3(ind(:,1))) );
  k3 = - M \ (nu*k*Un3 + Cj);

  Un4 = Un + deltat*k3;
% nonlinear terms
  Cj  = sixth * ( Un4(ind(:,3)) .* ( Un4(ind(:,2)) + Un4(ind(:,3)) ) - Un4(ind(:,1)) .* ( Un4(ind(:,2)) + Un4(ind(:,1)) ) );
% subgrid terms
  Cj = Cj - constant_sub * ( abs(Un4(ind(:,3))-Un4(ind(:,2))).*(Un4(ind(:,3))-Un4(ind(:,2))) - abs(Un4(ind(:,2))-Un4(ind(:,1))).*(Un4(ind(:,2))-Un4(ind(:,1))) );
  k4 = - M \ (nu*k*Un4 + Cj);

  y = Un + deltat*sixth*(k1 + 2*k2 + 2*k3 +k4 ) + F;
end
