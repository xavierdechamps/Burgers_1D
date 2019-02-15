function FE_HermiteH3(N,nu,constant_sub,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
% Solve the 1D forced Burgers equation with cubic Hermite elements 
% The unknowns are the velocity and the first spatial derivative of the velocity, thus 2 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element cubic Hermite H3")
  disp("************************************************************")

  h  = L/N; % Length of the elements
  DG = 4;   % Number of Gauss points for numerical integration of the energy
  X  = zeros(N,1);
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;

% ******** Mass Mij and rigidity Kij matrices *************
  Mij=h/35*[  13     11*h/6    9/2   -13*h/12  ;
             11*h/6   h*h/3   13*h/12  -h*h/4    ;
             9/2      13*h/12   13    -11*h/6;
            -13*h/12 -h*h/4   -11*h/6   h*h/3];
% LUMPED MATRIX
%  Mij = h * [  0.5     0       0      0     ;
%              0      0         0      0     ;
%              0      0        0.5     0     ;
%              0      0         0      0      ];

  Kij=[6/5/h 1/10 -6/5/h 1/10;
    1/10 2*h/15 -1/10 -h/30
    -6/5/h -1/10 6/5/h -1/10
    1/10 -h/30 -1/10 2*h/15];

	ind = zeros(N,6);
	ind(:,1) = -1:2:2*(N-1)-1; 
	ind(:,2) = 0:2:2*(N-1); 
	ind(:,3) = 1:2:2*(N-1)+1;
	ind(:,4) = 2:2:2*(N-1)+2;
	ind(:,5) = 3:2:2*(N-1)+3;
	ind(:,6) = 4:2:2*(N-1)+4;
	ind(1,:) = [2*N-1 2*N 1 2 3 4]; % Periodic condition
	ind(N,:) = [2*N-3 2*N-2 2*N-1 2*N 1 2]; % Periodic condition
  
% ************* Assemble the matrices (with periodic condition) ***********
  M=sparse(2*N,2*N);
  K=sparse(2*N,2*N);
  for i=1:2:2*N-2
    M(i:i+3,i:i+3)=M(i:i+3,i:i+3)+Mij;
    K(i:i+3,i:i+3)=K(i:i+3,i:i+3)+Kij;
  end
  M(2*N-1:2*N,2*N-1:2*N) = M(2*N-1:2*N,2*N-1:2*N) + Mij(1:2,1:2);
  M(2*N-1:2*N,1:2)       = M(2*N-1:2*N,1:2)       + Mij(1:2,3:4);
  M(1:2,1:2)             = M(1:2,1:2)             + Mij(3:4,3:4);
  M(1:2,2*N-1:2*N)       = M(1:2,2*N-1:2*N)       + Mij(3:4,1:2);

  K(2*N-1:2*N,2*N-1:2*N) = K(2*N-1:2*N,2*N-1:2*N) + Kij(1:2,1:2); 
  K(2*N-1:2*N,1:2)       = K(2*N-1:2*N,1:2)       + Kij(1:2,3:4);
  K(1:2,1:2)             = K(1:2,1:2)             + Kij(3:4,3:4);
  K(1:2,2*N-1:2*N)       = K(1:2,2*N-1:2*N)       + Kij(3:4,1:2);

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1)=2*rand(2*N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(1:2:2*N,1)=sin(X(:) * 2*pi/L) ; u(2:2:2*N,1)=2*pi/L*cos(X(:) * 2*pi/L) ;
  
  timeBeforeStatistics = 10;
  
  j=2;  z=2;  ind_error=1;
  filename=[name,num2str(1),'.mat'];
  uu=u(:,1);
%save(filename,'uu');

  kinEnergy = zeros(1:nbrpointtime+1,1); kinEnergy(1) = get_kinematic_energy(h,DG,u(:,1),N,2); kinEnergyMean = 0; nbrPointsStatistics = 0;
  spectralEnergy=zeros(Ninterpolation*N,1);
  reference_spectrum=load(file_spectrum);

  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2 = 2*pi*rand();    phi3 = 2*pi*rand();
    
    KForce = 2*pi/L;
    F(1:2:2*N,1)      =   sqrt(deltat) *         (   cos(2*KForce*X+phi2) +   cos(3*KForce*X+phi3) );
    F(2:2:2*N,1)      = - sqrt(deltat) * KForce* ( 2*sin(2*KForce*X+phi2) + 3*sin(3*KForce*X+phi3) );
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;    
  
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z) = RK4_FE_HermiteH3 (u(:,z-1),deltat,h,N,M,K,nu,F,constant_sub,ind);
    
    kinEnergy(i) = get_kinematic_energy(h,DG,u(:,end),N,2);
    
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1) ;
      Uinterpolated        = Interpolation(u(:,end),h,N,1,Ninterpolation); % Interpolate within each element to show internal oscillations
      length_Uinterpolated = length(Uinterpolated);
      fft_u                = fft(Uinterpolated);
      newSpectralEnergy    = fft_u.*conj(fft_u)/length_Uinterpolated/length_Uinterpolated;            
      spectralEnergy       = spectralEnergy + 2*pi*newSpectralEnergy;
      nbrPointsStatistics  = nbrPointsStatistics + 1;
    end
% Save the results, free some memory and show the results
    if ( mod(z/nbrpointtime,0.1) == 0)
        uu=u(:,end);
%        filename=[name,num2str(j),'.mat']; save(filename,'uu'); j=j+1;
        clear u;
        z=1;
        u(:,1)=uu;
        
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=['Spectral_energy_',name,'.mat']; spectralEnergyOut = spectralEnergy(1:N)/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
        
        Uinterpolated = Interpolation(u(:,end),h,N,1,Ninterpolation); % Interpolate within each element to show internal oscillations
        
        subplot(2,2,1)
        plot(0:1/(length(Uinterpolated)):1,[Uinterpolated;Uinterpolated(1)],'b','Linewidth',3)
        title(strcat('Time= ',num2str(i*deltat),', Re= ',num2str(mean(Uinterpolated)*L/nu)))
        xlabel('x/(2*\pi)'); ylabel('u(t)'); grid on
        xlim([0 1])
          
%          hold on
%          sol_theory = get_analytical_solution((i-1)*deltat,nu,[0:(2/(N*Ninterpolation)):1],100)' ;
%          plot([0:(2/(N*Ninterpolation)):1]/2,sol_theory,'r')
%          hold off
%          relative_error(ind_error,1) = (i-1)*deltat;
%          relative_error(ind_error,2) = sqrt( sum( (sol_theory-Uinterpolated(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%          disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%          ind_error = ind_error + 1;
          
%          hold on
%          plot([X;X(end)+h]/L,[u(1:2:2*N-1,end);u(1,end)],'b*')
%          hold off
        
        subplot(2,2,3);
        plot((1:i)*deltat,kinEnergy(1:i),'b','Linewidth',3);
        grid on; xlabel('Time'); ylabel('E(t)');
        xlim([0 time])
        
        subplot(1,2,2) ;
        loglog(0:N-1,spectralEnergy(1:N)/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
        grid on; xlabel('k'); ylabel('E(k)');
        xlim([1 reference_spectrum(end,1)])

        drawnow ;
    end
    
    z=z+1;
    CFL=u(1:2:end,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        break;
    end
  end

%  relative_error
  
  spectralEnergyOut = spectralEnergy(1:N)/nbrPointsStatistics;
  filename2=['Spectral_energy_',name,'.mat'];
  save(filename2,'spectralEnergyOut');

%  Uinterpolated = Interpolation(u(:,end),h,N,1,Ninterpolation);
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');

%  filename=strcat('Energy_',name,'.mat');
%  save(filename,'kinEnergy');
  
end

function y = RK4_FE_HermiteH3 (u,deltat,h,N,M,k,nu,F,constant_sub,ind)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
	
	Un = u ;
		
%	w   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % weight for numerical integration
%    ksi = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % coordinate of point for numerical integration
%    for i=1:length(w)
%		d_shape_fct_vector(:,i) = get_deriv_shape_fct(ksi(i));		
%	end
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Cj = get_non_linear_hermite_H3(Un,ind,N,h) ;
%	Cj += constant_sub * h *h * get_subgrid_terms(h,Un,ind,N,d_shape_fct_vector);

	k1 = - M \ (nu*k*Un + Cj);
	Un2 = Un + deltat*0.5*k1;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Cj = get_non_linear_hermite_H3(Un2,ind,N,h) ;
%	Cj += constant_sub * h *h * get_subgrid_terms(h,Un2,ind,N,d_shape_fct_vector);
	
	k2 = - M \ (nu*k*Un2 + Cj);
	Un3 = Un + deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  Cj = get_non_linear_hermite_H3(Un3,ind,N,h) ;
%	Cj += + constant_sub * h *h * get_subgrid_terms(h,Un3,ind,N,d_shape_fct_vector);
	
	k3 = - M \ (nu*k*Un3 + Cj);
	Un4 = Un + deltat*k3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  Cj = get_non_linear_hermite_H3(Un4,ind,N,h) ;
%	Cj += constant_sub * h *h * get_subgrid_terms(h,Un4,ind,N,d_shape_fct_vector);
	
	k4 = - M \ (nu*k*Un4 + Cj);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 + F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cnonlinear = get_non_linear_hermite_H3(Un,ind,N,h)
% Compute the discretization of the nonlinear convection term u * (du/dx)
	Cnonlinear = zeros(2*N,1);
	
	u1 = Un(ind(:,1));	du1 = Un(ind(:,2));  	u2 = Un(ind(:,3));  	du2 = Un(ind(:,4));  	u3 = Un(ind(:,5));  	du3 = Un(ind(:,6));
	
	Cnonlinear(1:2:2*N) = (du3.^2 - du1.^2)*h^2/168 + (du1.*du2 - du2.*du3)*h^2/105 - ...
                        5*(du1.*u1 + du3.*u3)*h/84 + (17*du2.*u1 - 17*du1.*u2)*h/420 + ...
                        (5*du2.*u2)*h/42 + 17*(du2.*u3 - du3.*u2)*h/420 + (u2.*u3 + u3.^2 -u1.^2 - u1.*u2)/6 ;

	Cnonlinear(2:2:2*N) = ( du1.^2 - du1.*du2 - du2.*du3 + du3.^2 )*h^3/840 + 11*h^2/840*(du1.*u1 - du3.*u3) + ...
                        17*(u1.^2 + u3.^2)*h/420 + (du1.*u2 - du3.*u2)*h^2/280 + 2*h*(u1.*u2 + u2.*u3)/105 - ...
                        (5*h*u2.^2)/42 + (du2.*u3 - du2.*u1)*h^2/168 ;
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
