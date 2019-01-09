function y = RK4_FE_HermiteP5 (u,deltat,h,N,M,k,nu,F,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with 5th order Hermite elements
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
  
  Un=u;

  ind = zeros(N,9);
  ind(:,1) = -2:3:3*N-5; 
  ind(:,2) = -1:3:3*N-4; 
  ind(:,3) = 0:3:3*N-3;
  ind(:,4) = 1:3:3*N-2;
  ind(:,5) = 2:3:3*N-1;
  ind(:,6) = 3:3:3*N;
  ind(:,7) = 4:3:3*N+1;
  ind(:,8) = 5:3:3*N+2;
  ind(:,9) = 6:3:3*N+3;
  ind(1,:) = [3*N-2 3*N-1 3*N 1 2 3 4 5 6]; 
  ind(N,:) = [3*N-5 3*N-4 3*N-3 3*N-2 3*N-1 3*N 1 2 3];
  
%  w   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % weight for numerical integration
%    ksi = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % coordinate of point for numerical integration
%    for i=1:length(w)
%    d_shape_fct_vector(:,i) = get_deriv_shape_fct(ksi(i));    
%  end

  Cj = get_non_linear_hermite_P5(Un,ind,N,h) ;
    
  k1 = - M \ (nu*k*Un + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un2 = Un + deltat*0.5*k1;
  
  Cj = get_non_linear_hermite_P5(Un2,ind,N,h) ;
    
%  Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un2,ind,N,d_shape_fct_vector);
  
  k2 = - M \ (nu*k*Un2 + Cj);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un3 = Un + deltat*0.5*k2;
  
  Cj = get_non_linear_hermite_P5(Un3,ind,N,h) ;
    
%  Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un3,ind,N,d_shape_fct_vector);
  
  k3 = - M \ (nu*k*Un3 + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un4 = Un + deltat*k3;

  Cj = get_non_linear_hermite_P5(Un4,ind,N,h) ;
    
%  Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un4,ind,N,d_shape_fct_vector);
  
  k4 = - M \ (nu*k*Un4 + Cj);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y = Un + deltat/6*( k1 + 2*k2 + 2*k3 + k4 ) + F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cnonlinear = get_non_linear_hermite_P5(Un,ind,N,h)
% Compute the discretization of the nonlinear convection term u * (du/dx)
  Cnonlinear = zeros(3*N,1);
    
  u1 = Un(ind(:,1));  du1 = Un(ind(:,2));    ddu1 = Un(ind(:,3));
  u2 = Un(ind(:,4));  du2 = Un(ind(:,5));    ddu2 = Un(ind(:,6));
  u3 = Un(ind(:,7));  du3 = Un(ind(:,8));    ddu3 = Un(ind(:,9));
  
  Cnonlinear(1:3:3*N) = (-496*h^2*du1.^2 + 732*du1.*du2*h^2 - 732*du2.*du3*h^2 + 496*h^2*du3.^2 - 88*ddu1.*du1*h^3 - 72*ddu2.*du1*h^3 + 72*ddu1.*du2*h^3 + 176*ddu2.*du2*h^3 + 72*ddu3.*du2*h^3 - 72*ddu2.*du3*h^3 - 88*ddu3.*du3*h^3 - 4*h^4*ddu1.^2 - 7*ddu1.*ddu2*h^4 + 7*ddu2.*ddu3*h^4 + 4*h^4*ddu3.^2 - 3848*h*du1.*u1 + 2444*h*du2.*u1 - 328*h^2*ddu1.*u1 - 244*h^2*ddu2.*u1 - 8008*u1.^2 - 2444*h*du1.*u2 + 7696*h*du2.*u2 - 2444*h*du3.*u2 - 244*h^2*ddu1.*u2 + 244*h^2*ddu3.*u2 - 8008*u1.*u2 + 2444*h*du2.*u3 - 3848*h*du3.*u3 + 244*h^2*ddu2.*u3 + 328*h^2*ddu3.*u3 + 8008*u2.*u3 + 8008*u3.^2) /48048 ;
  Cnonlinear(2:3:3*N) = (2032*h^3*du1.^2 - 2032*du1.*du2*h^3 - 2032*du2.*du3*h^3 + 2032*h^3*du3.^2 + 345*ddu1.*du1*h^4 + 220*ddu2.*du1*h^4 - 184*ddu1.*du2*h^4 + 184*ddu3.*du2*h^4 - 220*ddu2.*du3*h^4 - 345*ddu3.*du3*h^4 + 15*ddu1.^2*h^5 + 20*ddu1.*ddu2*h^5 + 12*ddu2.^2*h^5 + 20*ddu2.*ddu3*h^5 + 15*ddu3.^2*h^5 + 16644*h^2*du1.*u1 - 7440*h^2*du2.*u1 + 1357*h^3*ddu1.*u1 + 805*h^3*ddu2.*u1 + 36660*h*u1.^2 + 5664*h^2*du1.*u2 - 5664*h^2*du3.*u2 + 502*h^3*ddu1.*u2 - 180*h^3*ddu2.*u2 + 502*h^3*ddu3.*u2 + 21060*h*u1.*u2 - 115440*h*u2.^2 + 7440*h^2*du2.*u3 - 16644*h^2*du3.*u3 + 805*h^3*ddu2.*u3 + 1357*h^3*ddu3.*u3 + 21060*h*u2.*u3 + 36660*h*u3.^2) /720720 ;
  Cnonlinear(3:3:3*N) = (-736*du1.^2*h^4 + 500*du1.*du2*h^4 - 500*du2.*du3*h^4 + 736*du3.^2*h^4 - 120*ddu1.*du1*h^5 - 60*ddu2.*du1*h^5 + 40*ddu1.*du2*h^5 - 48*ddu2.*du2*h^5 + 40*ddu3.*du2*h^5 - 60*ddu2.*du3*h^5 - 120*ddu3.*du3*h^5 - 5*ddu1.^2*h^6 - 5*ddu1.*ddu2*h^6 + 5*ddu2.*ddu3*h^6 + 5*ddu3.^2*h^6 - 6328*h^3*du1.*u1 + 2060*h^3*du2.*u1 - 496*h^4*ddu1.*u1 - 240*h^4*ddu2.*u1 - 14640*h^2*u1.^2 - 1108*h^3*du1.*u2 - 9840*h^3*du2.*u2 - 1108*h^3*du3.*u2 - 76*h^4*ddu1.*u2 + 76*h^4*ddu3.*u2 - 5040*h^2*u1.*u2 + 2060*h^3*du2.*u3 - 6328*h^3*du3.*u3 + 240*h^4*ddu2.*u3 + 496*h^4*ddu3.*u3 +  5040*h^2*u2.*u3 + 14640*h^2*u3.^2) /2882880 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%