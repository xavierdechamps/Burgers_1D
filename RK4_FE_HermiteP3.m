function y = RK4_FE_HermiteP3 (u,deltat,h,N,M,k,nu,F,constant_sub)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with cubic Hermite elements
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

	ind = zeros(N,6);
	ind(:,1) = -1:2:2*(N-1)-1; 
	ind(:,2) = 0:2:2*(N-1); 
	ind(:,3) = 1:2:2*(N-1)+1;
	ind(:,4) = 2:2:2*(N-1)+2;
	ind(:,5) = 3:2:2*(N-1)+3;
	ind(:,6) = 4:2:2*(N-1)+4;
	ind(1,:) = [2*N-1 2*N 1 2 3 4]; % Periodic condition
	ind(N,:) = [2*N-3 2*N-2 2*N-1 2*N 1 2]; % Periodic condition
		
%	w   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % weight for numerical integration
%    ksi = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % coordinate of point for numerical integration
%    for i=1:length(w)
%		d_shape_fct_vector(:,i) = get_deriv_shape_fct(ksi(i));		
%	end
			
  Cj = get_non_linear_hermite_P3(Un,ind,N,h) ;
    
%	Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un,ind,N,d_shape_fct_vector);

	k1 = - M \ (nu*k*Un + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Un2 = Un + deltat*0.5*k1;
	
  Cj = get_non_linear_hermite_P3(Un2,ind,N,h) ;
    
%	Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un2,ind,N,d_shape_fct_vector);
	
	k2 = - M \ (nu*k*Un2 + Cj);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Un3 = Un + deltat*0.5*k2;
	
  Cj = get_non_linear_hermite_P3(Un3,ind,N,h) ;
    
%	Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un3,ind,N,d_shape_fct_vector);
	
	k3 = - M \ (nu*k*Un3 + Cj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Un4 = Un + deltat*k3;

  Cj = get_non_linear_hermite_P3(Un4,ind,N,h) ;
    
%	Cj = Cj + constant_sub * h *h * get_subgrid_terms(h,Un4,ind,N,d_shape_fct_vector);
	
	k4 = - M \ (nu*k*Un4 + Cj);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y = Un + deltat/6*( k1 + 2*k2 + 2*k3 + k4 ) + F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cnonlinear = get_non_linear_hermite_P3(Un,ind,N,h)
% Compute the discretization of the nonlinear convection term u * (du/dx)
	Cnonlinear = zeros(2*N,1);
	
	u1 = Un(ind(:,1));	du1 = Un(ind(:,2));  	u2 = Un(ind(:,3));  	du2 = Un(ind(:,4));  	u3 = Un(ind(:,5));  	du3 = Un(ind(:,6));
	
	Cnonlinear(1:2:2*N) = (du3.^2 - du1.^2)*h^2/168 + (du1.*du2 - du2.*du3)*h^2/105 - 5*(du1.*u1 + du3.*u3)*h/84 + (17*du2.*u1 - 17*du1.*u2)*h/420 + (5*du2.*u2)*h/42 + 17*(du2.*u3 - du3.*u2)*h/420 + (u2.*u3 + u3.^2 -u1.^2 - u1.*u2)/6 ;

	Cnonlinear(2:2:2*N) = ( du1.^2 - du1.*du2 - du2.*du3 + du3.^2 )*h^3/840 + 11*h^2/840*(du1.*u1 - du3.*u3) + 17*(u1.^2 + u3.^2)*h/420 + (du1.*u2 - du3.*u2)*h^2/280 + 2*h*(u1.*u2 + u2.*u3)/105 - (5*h*u2.^2)/42 + (du2.*u3 - du2.*u1)*h^2/168 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d_shape_fct = get_deriv_shape_fct(ksi)
% Analytical expression of the derivative of the shape functions
	d_shape_fct = [ 0.5 * (ksi - 2)*(1 + ksi) + 0.25 * (1 + ksi)^2                             ;
	                0.5 * (1 + ksi)*((1 + ksi)*0.5 - 1) + 0.5* (0.25* (1 + ksi)^2 -ksi)        ;
	                0.5 * (2 - ksi)*(1 + ksi) - 0.25* (1 + ksi)^2                              ;
	                0.25* (ksi -1)*(1 + ksi) + 0.125* (1 + ksi)^2        
	              ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Csubgrid = get_subgrid_terms(h,Un,ind,N,deriv_shape)
% Compute the subgrid terms needed to model turbulence
	Csubgrid = zeros(2*N,1);
	
	u1 = Un(ind(:,1));	du1 = Un(ind(:,2));  	u2 = Un(ind(:,3));  	du2 = Un(ind(:,4));  	u3 = Un(ind(:,5));  	du3 = Un(ind(:,6));

	Csubgrid(1:2:2*N) = 103*(du3.^2 - du1.^2)/8505 + 113*(du2.*du3 - du1.*du2)/2835 - 3904*(du1.*u1 + du2.*u1 - du3.*u2)/(8505*h) - 54*(u1.^2 - 2*u2.^2 + u3.^2)/(35*h^2) + 988*(du1.*u2 - du2.*u3 - du3.*u3)/(8505*h) + 4892*du2.*u2/(8505*h) ;

	Csubgrid(2:2:2*N) = -113*(du1.^2 + du3.^2)*h/8505 + 452*(du2.^2)*h/2835 + 2*(du1.*du2 + du2.*du3)*h/315 + 16*(du1.*u1 + du3.*u2)/2835 + (832*du2.*u1)/8505 + 6*(u1.^2 - u3.^2)/(35*h) - 97*(du1.*u2 + du3.*u3)/2835 + (6*du2.*u2)/35 + (626*du2.*u3)/8505 ;

	
%	for i = 1:size(deriv_shape,2)
	
%		deriv_u = deriv_shape(1,i) * u1 + deriv_shape(2,i) * du1 + deriv_shape(3,i) * u2 + deriv_shape(4,i) * du2 ;
				
%		Csubgrid(ind(:,3)) = Csubgrid(ind(:,3)) + 0.125*deriv_shape(1,i)*deriv_u.*abs(deriv_u);
%		Csubgrid(ind(:,4)) = Csubgrid(ind(:,4)) + 0.125*deriv_shape(2,i)*deriv_u.*abs(deriv_u);
%		Csubgrid(ind(:,5)) = Csubgrid(ind(:,5)) + 0.125*deriv_shape(3,i)*deriv_u.*abs(deriv_u);
%		Csubgrid(ind(:,6)) = Csubgrid(ind(:,6)) + 0.125*deriv_shape(4,i)*deriv_u.*abs(deriv_u);
%	end
	
end
