function y = RK4_FE_Lagrangep3 (u,deltat,h,N,M,k,nu,F,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with cubic Lagrange elements
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

  Cj=zeros(3*N,1);
  sixth = 1/6;
  half = 0.5;

  ind = zeros(N,4);
  ind(:,1) = 1:3:3*N-2; 
  ind(:,2) = 2:3:3*N-1; 
  ind(:,3) = 3:3:3*N;
  ind(:,4) = 4:3:3*N+1;
  ind(end,4) = 1;

  w   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % weight for numerical integration
  ksi = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % coordinate of point for numerical integration
  for i=1:length(w)
    d_shape_fct_vector(:,i) = get_deriv_shape_fct(ksi(i));    
  end
  
  Un = u;
  u1 = Un(ind(:,1));  u2 = Un(ind(:,2));    u3 = Un(ind(:,3));    u4 = Un(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (constant_sub~=0)
    Cj = Cj + get_subgrid_terms(constant_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector) ;
  end
% RK step 1
  k1 = - M \ (nu*k*Un + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un2 = Un + deltat*half*k1;
  u1  = Un2(ind(:,1));  u2 = Un2(ind(:,2));    u3 = Un2(ind(:,3));    u4 = Un2(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (constant_sub~=0)
    Cj = Cj + get_subgrid_terms(constant_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector) ;
  end
% RK step 2
  k2 = - M \ (nu*k*Un2 + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un3 = Un + deltat*half*k2;
  u1  = Un3(ind(:,1));  u2 = Un3(ind(:,2));    u3 = Un3(ind(:,3));    u4 = Un3(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (constant_sub~=0)
    Cj = Cj + get_subgrid_terms(constant_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector) ;
  end
% RK step 3
  k3 = - M \ (nu*k*Un3 + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Un4 = Un + deltat*k3;
  u1  = Un4(ind(:,1));  u2 = Un4(ind(:,2));    u3 = Un4(ind(:,3));    u4 = Un4(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (constant_sub~=0)
    Cj = Cj + get_subgrid_terms(constant_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector) ;
  end
% RK step 4
  k4 = - M \ (nu*k*Un4 + Cj);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK step 5
  y = Un + deltat*sixth*(k1 + 2*k2 + 2*k3 +k4 ) + F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cnonlinear = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N)
% Compute the discretization of the nonlinear convection term u * (du/dx)
  Cnonlinear = zeros(3*N,1);
  
  Cnonlinear(ind(:,1)) = Cnonlinear(ind(:,1)) + (-(u1.^2/3) + (537*u1.*u2)/2240 + (135*u2.^2)/448 - (3*u1.*u3)/32 - ( 351*u2.*u3)/2240 - (27*u3.^2)/1120 +(139*u1.*u4)/6720 + ( 3*u2.*u4)/112 - (3*u3.*u4)/2240 + (139*u4.^2)/6720);
  Cnonlinear(ind(:,2)) = Cnonlinear(ind(:,2)) + (-((537*u1.^2)/2240) - (135*u1.*u2)/448 + (27*u1.*u3)/280 + ( 729*u2.*u3)/2240 + (729*u3.^2)/2240 - (9*u1.*u4)/320 - ( 27*u2.*u4)/1120 - (27*u3.*u4)/448 - (3*u4.^2)/32);
  Cnonlinear(ind(:,3)) = Cnonlinear(ind(:,3)) + ((3*u1.^2)/32 + (27*u1.*u2)/448 - (729*u2.^2)/2240 + (27*u1.*u3)/1120 - ( 729*u2.*u3)/2240 + (9*u1.*u4)/320 - (27*u2.*u4)/280 + ( 135*u3.*u4)/448 + (537*u4.^2)/2240);
  Cnonlinear(ind(:,4)) = Cnonlinear(ind(:,4)) + (-((139*u1.^2)/6720) + (3*u1.*u2)/2240 + (27*u2.^2)/1120 - ( 3*u1.*u3)/112 + (351*u2.*u3)/2240 - (135*u3.^2)/448 - ( 139*u1.*u4)/6720 + (3*u2.*u4)/32 - (537*u3.*u4)/2240 + u4.^2/3);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_shape_fct = get_deriv_shape_fct(ksi)
% Analytical expression of the derivative of the shape functions
  d_shape_fct = [ (1 - ksi) * ksi * 9/8      + (1 - 9*ksi^2)/16            ;
                  ksi * (3 * ksi -1) * 9/8   + (ksi^2 - 1) * 27/16         ;
                  (1 - ksi^2) * 27/16        - ksi * (1 + 3 * ksi) * 9/8   ;
                  ksi * (1 + ksi) * 9/8      + (9 * ksi^2 - 1) / 16        ;
                ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Csubgrid = get_subgrid_terms(C,h,u1,u2,u3,u4,ind,N,deriv_shape)
% Compute the subgrid terms needed to model turbulence
  Csubgrid = zeros(3*N,1);
  eightth = 0.125;
  for i = 1:size(deriv_shape,2)
    deriv_u = deriv_shape(1,i) * u1 + deriv_shape(2,i) * u2 + deriv_shape(3,i) * u3 + deriv_shape(4,i) * u4 ;
%    Csubgrid(ind(:,1)) = Csubgrid(ind(:,1)) + C*h*h*h*0.5*deriv_shape(1,i)*deriv_u.*abs(deriv_u);%
%    Csubgrid(ind(:,2)) = Csubgrid(ind(:,2)) + C*h*h*h*0.5*deriv_shape(2,i)*deriv_u.*abs(deriv_u);
%    Csubgrid(ind(:,3)) = Csubgrid(ind(:,3)) + C*h*h*h*0.5*deriv_shape(3,i)*deriv_u.*abs(deriv_u);
%    Csubgrid(ind(:,4)) = Csubgrid(ind(:,4)) + C*h*h*h*0.5*deriv_shape(4,i)*deriv_u.*abs(deriv_u);
    
    Csubgrid(ind(:,1)) = Csubgrid(ind(:,1)) + C*eightth*deriv_shape(1,i)*deriv_u.*abs(deriv_u);
    Csubgrid(ind(:,2)) = Csubgrid(ind(:,2)) + C*eightth*deriv_shape(2,i)*deriv_u.*abs(deriv_u);
    Csubgrid(ind(:,3)) = Csubgrid(ind(:,3)) + C*eightth*deriv_shape(3,i)*deriv_u.*abs(deriv_u);
    Csubgrid(ind(:,4)) = Csubgrid(ind(:,4)) + C*eightth*deriv_shape(4,i)*deriv_u.*abs(deriv_u);
  end
  
end
