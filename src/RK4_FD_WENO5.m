function y = RK4_FD_WENO5 (u,deltat,N,K,F,h,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with a nonlinear discretization WENO for the convective term
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
  ind = zeros(N,7);
  ind(:,1) = -2:N-3; 
  ind(:,2) = -1:N-2; 
  ind(:,3) = 0:N-1;
  ind(:,4) = 1:N;
  ind(:,5) = 2:N+1;
  ind(:,6) = 3:N+2;
  ind(:,7) = 4:N+3;
  ind(1,:) = [N-2 N-1 N 1 2 3 4]; 
  ind(2,:) = [N-1  N  1 2 3 4 5]; 
  ind(3,:) = [ N   1  2 3 4 5 6]; 
  ind(N-2,:) = [N-5 N-4 N-3 N-2 N-1 N 1];
  ind(N-1,:) = [N-4 N-3 N-2 N-1  N  1 2];
  ind(N  ,:) = [N-3 N-2 N-1  N   1  2 3];

  Un = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step  
  C   = get_non_linear_WENO5(Un,ind,N,h) ;
  k1  = K * Un + C;
  Un2 = Un + deltat * 0.5 * k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C   = get_non_linear_WENO5(Un2,ind,N,h) ;
  k2  = K * Un2 + C;
  Un3 = Un + deltat * 0.5 * k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step  
  C   = get_non_linear_WENO5(Un3,ind,N,h) ;
  k3  = K * Un3 + C;
  Un4 = Un + deltat * k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step  
  C   = get_non_linear_WENO5(Un4,ind,N,h) ;
  k4 = K * Un4 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y  = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 +F;
  
end

%---------------------------------------------------------------------------------------------------------------
function C   = get_non_linear_WENO5(u,i,N,h)
% i = i-3  i-2  i-1  i  i+1  i+2  i+3
%      1    2    3   4   5    6    7
  C = zeros(N,1);
  epsi = 1e-6;
  
  for j=1:7
    matu(:,j) = u(i(:,j));
  end
  theta = max(matu,[],2);
% beta for node i+1/2
  betapip1(1:N,1) = (13/12)*( u(i(:,2)).*(0.5*u(i(:,2))+theta) - 2*u(i(:,3)).*(0.5*u(i(:,3))+theta) + u(i(:,4)).*(0.5*u(i(:,4))+theta)).^2 + 0.25*( u(i(:,2)).*(0.5*u(i(:,2))+theta) - 4*u(i(:,3)).*(0.5*u(i(:,3))+theta) + 3*u(i(:,4)).*(0.5*u(i(:,4))+theta) ).^2 ;
  betapip1(1:N,2) = (13/12)*( u(i(:,3)).*(0.5*u(i(:,3))+theta) - 2*u(i(:,4)).*(0.5*u(i(:,4))+theta) + u(i(:,5)).*(0.5*u(i(:,5))+theta)).^2 + 0.25*( u(i(:,3)).*(0.5*u(i(:,3))+theta) - u(i(:,5)).*(0.5*u(i(:,5))+theta) ).^2 ;
  betapip1(1:N,3) = (13/12)*( u(i(:,4)).*(0.5*u(i(:,4))+theta) - 2*u(i(:,5)).*(0.5*u(i(:,5))+theta) + u(i(:,6)).*(0.5*u(i(:,6))+theta)).^2 + 0.25*( 3*u(i(:,4)).*(0.5*u(i(:,4))+theta) - 4*u(i(:,5)).*(0.5*u(i(:,5))+theta) + u(i(:,6)).*(0.5*u(i(:,6))+theta) ).^2 ;
  betamip1(1:N,1) = (13/12)*( u(i(:,5)).*(0.5*u(i(:,5))-theta) - 2*u(i(:,6)).*(0.5*u(i(:,6))-theta) + u(i(:,7)).*(0.5*u(i(:,7))-theta)).^2 + 0.25*( 3*u(i(:,5)).*(0.5*u(i(:,5))-theta) - 4*u(i(:,6)).*(0.5*u(i(:,6))-theta) + u(i(:,7)).*(0.5*u(i(:,7))-theta) ).^2 ;
  betamip1(1:N,2) = (13/12)*( u(i(:,4)).*(0.5*u(i(:,4))-theta) - 2*u(i(:,5)).*(0.5*u(i(:,5))-theta) + u(i(:,6)).*(0.5*u(i(:,6))-theta)).^2 + 0.25*( u(i(:,4)).*(0.5*u(i(:,4))-theta) - u(i(:,6)).*(0.5*u(i(:,6))-theta) ).^2 ;
  betamip1(1:N,3) = (13/12)*( u(i(:,3)).*(0.5*u(i(:,3))-theta) - 2*u(i(:,4)).*(0.5*u(i(:,4))-theta) + u(i(:,5)).*(0.5*u(i(:,5))-theta)).^2 + 0.25*( u(i(:,3)).*(0.5*u(i(:,3))-theta) - 4*u(i(:,4)).*(0.5*u(i(:,4))-theta) + 3*u(i(:,5)).*(0.5*u(i(:,5))-theta) ).^2 ;
% beta for node i-1/2
  betapim1(1:N,1) = (13/12)*( u(i(:,1)).*(0.5*u(i(:,1))+theta) - 2*u(i(:,2)).*(0.5*u(i(:,2))+theta) + u(i(:,3)).*(0.5*u(i(:,3))+theta)).^2 + 0.25*( u(i(:,1)).*(0.5*u(i(:,1))+theta) - 4*u(i(:,2)).*(0.5*u(i(:,2))+theta) + 3*u(i(:,3)).*(0.5*u(i(:,3))+theta) ).^2 ;
  betapim1(1:N,2) = (13/12)*( u(i(:,2)).*(0.5*u(i(:,2))+theta) - 2*u(i(:,3)).*(0.5*u(i(:,3))+theta) + u(i(:,4)).*(0.5*u(i(:,4))+theta)).^2 + 0.25*( u(i(:,2)).*(0.5*u(i(:,2))+theta) - u(i(:,4)).*(0.5*u(i(:,4))+theta) ).^2 ;
  betapim1(1:N,3) = (13/12)*( u(i(:,3)).*(0.5*u(i(:,3))+theta) - 2*u(i(:,4)).*(0.5*u(i(:,4))+theta) + u(i(:,5)).*(0.5*u(i(:,5))+theta)).^2 + 0.25*( 3*u(i(:,3)).*(0.5*u(i(:,3))+theta) - 4*u(i(:,4)).*(0.5*u(i(:,4))+theta) + u(i(:,5)).*(0.5*u(i(:,5))+theta) ).^2 ;
  betamim1(1:N,1) = (13/12)*( u(i(:,4)).*(0.5*u(i(:,4))-theta) - 2*u(i(:,5)).*(0.5*u(i(:,5))-theta) + u(i(:,6)).*(0.5*u(i(:,6))-theta)).^2 + 0.25*( 3*u(i(:,4)).*(0.5*u(i(:,4))-theta) - 4*u(i(:,5)).*(0.5*u(i(:,5))-theta) + u(i(:,6)).*(0.5*u(i(:,6))-theta) ).^2 ;
  betamim1(1:N,2) = (13/12)*( u(i(:,3)).*(0.5*u(i(:,3))-theta) - 2*u(i(:,4)).*(0.5*u(i(:,4))-theta) + u(i(:,5)).*(0.5*u(i(:,5))-theta)).^2 + 0.25*( u(i(:,3)).*(0.5*u(i(:,3))-theta) - u(i(:,5)).*(0.5*u(i(:,5))-theta) ).^2 ;
  betamim1(1:N,3) = (13/12)*( u(i(:,2)).*(0.5*u(i(:,2))-theta) - 2*u(i(:,3)).*(0.5*u(i(:,3))-theta) + u(i(:,4)).*(0.5*u(i(:,4))-theta)).^2 + 0.25*( u(i(:,2)).*(0.5*u(i(:,2))-theta) - 4*u(i(:,3)).*(0.5*u(i(:,3))-theta) + 3*u(i(:,4)).*(0.5*u(i(:,4))-theta) ).^2 ;
  
% alpha
  alphapip1 = [0.1./(epsi+betapip1(:,1))  0.6./(epsi+betapip1(:,2)) 0.3./(epsi+betapip1(:,3)) ] ;
  alphamip1 = [0.1./(epsi+betamip1(:,1))  0.6./(epsi+betamip1(:,2)) 0.3./(epsi+betamip1(:,3)) ] ;
  alphapim1 = [0.1./(epsi+betapim1(:,1))  0.6./(epsi+betapim1(:,2)) 0.3./(epsi+betapim1(:,3)) ] ;
  alphamim1 = [0.1./(epsi+betamim1(:,1))  0.6./(epsi+betamim1(:,2)) 0.3./(epsi+betamim1(:,3)) ] ;
% weight function
  wpip1 = alphapip1./sum(alphapip1,2) ; % The sum over the rows is equal to 1 ---> correct
  wmip1 = alphamip1./sum(alphamip1,2) ;
  wpim1 = alphapim1./sum(alphapim1,2) ;
  wmim1 = alphamim1./sum(alphamim1,2) ;
% fluxes
  fpip1 = wpip1(:,1).*(2*u(i(:,2)).*(0.5*u(i(:,2))+theta) - 7*u(i(:,3)).*(0.5*u(i(:,3))+theta) + 11*u(i(:,4)).*(0.5*u(i(:,4))+theta) ) + wpip1(:,2).*(-u(i(:,3)).*(0.5*u(i(:,3))+theta) + 5*u(i(:,4)).*(0.5*u(i(:,4))+theta) + 2*u(i(:,5)).*(0.5*u(i(:,5))+theta)) + wpip1(:,3).*(2*u(i(:,4)).*(0.5*u(i(:,4))+theta)+5*u(i(:,5)).*(0.5*u(i(:,5))+theta)-u(i(:,6)).*(0.5*u(i(:,6))+theta));
  fpim1 = wpim1(:,1).*(2*u(i(:,1)).*(0.5*u(i(:,1))+theta) - 7*u(i(:,2)).*(0.5*u(i(:,2))+theta) + 11*u(i(:,3)).*(0.5*u(i(:,3))+theta) ) + wpim1(:,2).*(-u(i(:,2)).*(0.5*u(i(:,2))+theta) + 5*u(i(:,3)).*(0.5*u(i(:,3))+theta) + 2*u(i(:,4)).*(0.5*u(i(:,4))+theta)) + wpim1(:,3).*(2*u(i(:,3)).*(0.5*u(i(:,3))+theta)+5*u(i(:,4)).*(0.5*u(i(:,4))+theta)-u(i(:,5)).*(0.5*u(i(:,5))+theta));
  fmip1 = wmip1(:,1).*(11*u(i(:,5)).*(0.5*u(i(:,5))-theta) - 7*u(i(:,6)).*(0.5*u(i(:,6))-theta) + 2*u(i(:,7)).*(0.5*u(i(:,7))-theta) ) + wmip1(:,2).*(2*u(i(:,4)).*(0.5*u(i(:,4))-theta) + 5*u(i(:,5)).*(0.5*u(i(:,5))-theta) - u(i(:,6)).*(0.5*u(i(:,6))-theta)) + wmip1(:,3).*(-u(i(:,3)).*(0.5*u(i(:,3))-theta)+5*u(i(:,4)).*(0.5*u(i(:,4))-theta)+2*u(i(:,5)).*(0.5*u(i(:,5))-theta));
  fmim1 = wmim1(:,1).*(11*u(i(:,4)).*(0.5*u(i(:,4))-theta) - 7*u(i(:,5)).*(0.5*u(i(:,5))-theta) + 2*u(i(:,6)).*(0.5*u(i(:,6))-theta) ) + wmim1(:,2).*(2*u(i(:,3)).*(0.5*u(i(:,3))-theta) + 5*u(i(:,4)).*(0.5*u(i(:,4))-theta) - u(i(:,5)).*(0.5*u(i(:,5))-theta)) + wmim1(:,3).*(-u(i(:,2)).*(0.5*u(i(:,2))-theta)+5*u(i(:,3)).*(0.5*u(i(:,3))-theta)+2*u(i(:,4)).*(0.5*u(i(:,4))-theta));
  
  C = - (fpip1 - fpim1 + fmip1 - fmim1)/(12*h) ;
end
