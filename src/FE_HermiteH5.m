function FE_HermiteH5(N,nu,constant_sub,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
% Solve the 1D forced Burgers equation with 5th order Hermite elements 
% The unknowns are the velocity, the first and second spatial derivatives of the velocity,
% thus 3 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element 5th order Hermite H5")
  disp("************************************************************")

  length_vec = 3*N;
  h  = L/N; % Length of the elements
  DG = 4;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;

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
  
% ******** Mass Mij and rigidity Kij matrices *************
  Mij= [ 181*h/462   311*(h^2)/4620   281*(h^3)/55440   25*h/231   -151*(h^2)/4620   181*(h^3)/55440 ;
       311*(h^2)/4620  52*(h^3)/3465  23*(h^4)/18480  151*(h^2)/4620  -19*(h^3)/1980  13*(h^4)/13860 ;
       281*(h^3)/55440  23*(h^4)/18480  (h^5)/9240  181*(h^3)/55440  -13*(h^4)/13860  (h^5)/11088 ;
       25*h/231  151*(h^2)/4620  181*(h^3)/55440  181*h/462  -311*(h^2)/4620  281*(h^3)/55440 ;
       -151*(h^2)/4620  -19*(h^3)/1980  -13*(h^4)/13860  -311*(h^2)/4620  52*(h^3)/3465  -23*(h^4)/18480 ;
       181*(h^3)/55440  13*(h^4)/13860  (h^5)/11088  281*(h^3)/55440  -23*(h^4)/18480  (h^5)/9240
  ];
  
  Kij=[10/(7*h)  3/14  h/84  -10/(7*h)  3/14  -h/84  ;
     3/14  8*h/35  (h^2)/60  -3/14  -h/70  (h^2)/210 ;
     h/84  (h^2)/60  (h^3)/630  -h/84  -(h^2)/210  (h^3)/1260 ;
     -10/(7*h) -3/14  -h/84  10/(7*h) -3/14  h/84 ;
     3/14  -h/70  -(h^2)/210  -3/14  8*h/35  -(h^2)/60 ;
     -h/84  (h^2)/210  (h^3)/1260  h/84  -(h^2)/60  (h^3)/630
  ];

% ************* Assemble the matrices (with periodic condition) ***********
  M=sparse(length_vec,length_vec);
  K=sparse(length_vec,length_vec);
  for i=1:3:length_vec-3
    M(i:i+5,i:i+5)=M(i:i+5,i:i+5)+Mij;
    K(i:i+5,i:i+5)=K(i:i+5,i:i+5)+Kij;
  end
  M(3*N-2:3*N,3*N-2:3*N) = M(3*N-2:3*N,3*N-2:3*N) + Mij(1:3,1:3);
  M(3*N-2:3*N,1:3)       = M(3*N-2:3*N,1:3)       + Mij(1:3,4:6);
  M(1:3,1:3)             = M(1:3,1:3)             + Mij(4:6,4:6);
  M(1:3,3*N-2:3*N)       = M(1:3,3*N-2:3*N)       + Mij(4:6,1:3);

  K(3*N-2:3*N,3*N-2:3*N) = K(3*N-2:3*N,3*N-2:3*N) + Kij(1:3,1:3);
  K(3*N-2:3*N,1:3)       = K(3*N-2:3*N,1:3)       + Kij(1:3,4:6);
  K(1:3,1:3)             = K(1:3,1:3)             + Kij(4:6,4:6);
  K(1:3,3*N-2:3*N)       = K(1:3,3*N-2:3*N)       + Kij(4:6,1:3);

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1)=2*rand(length_vec,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(1:3:3*N,1)= sin(X(:) * 2*pi/L)              ;
%  u(2:3:3*N,1)= cos(X(:) * 2*pi/L) * (2*pi/L)   ;
%  u(3:3:3*N,1)=-sin(X(:) * 2*pi/L) * (2*pi/L)^2 ;

  timeBeforeStatistics = 10;
  
% kinematic energy not computed yet for the 5th order Hermite element
%  kinEnergy = zeros(1:nbrpointtime+1,1); kinEnergy(1)=get_kinematic_energy(h,DG,u(:,1),N,2);
  j=2;  z=2;  ind_error = 1;
  filename=[name,num2str(1),'.mat'];
  uu=u(:,1);
%  save(filename,'uu');

  nbrPointsStatistics=0;
  spectralEnergy=zeros(Ninterpolation*N,1);
  reference_spectrum=load(file_spectrum);

  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand();    phi3=2*pi*rand();
    KForce=2*pi/L;
    F(1:3:length_vec,1)      =   sqrt(deltat)                   * (   cos(2*KForce*X+phi2) +   cos(3*KForce*X+phi3) );
    F(2:3:length_vec,1)      = - sqrt(deltat) * KForce          * ( 2*sin(2*KForce*X+phi2) + 3*sin(3*KForce*X+phi3) );
    F(3:3:length_vec,1)      = - sqrt(deltat) * KForce * KForce * ( 4*cos(2*KForce*X+phi2) + 9*cos(3*KForce*X+phi3) );
% Uncomment the following line in case of sinus wave with non-linear convection
%      F = 0;    
  
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z) = RK4_FE_HermiteH5 (u(:,z-1),deltat,h,N,M,K,nu,F,constant_sub,ind);
% kinematic energy not computed yet for the 5th order Hermite element
%    kinEnergy(i)=CalculEnergie(h,DG,u(:,end),N,2);
    
    if (i*deltat>=timeBeforeStatistics)
      Uinterpolated        = Interpolation(u(:,end),h,N,3,Ninterpolation);% Interpolate within each element to show internal oscillations
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
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
      
        Uinterpolated = Interpolation(u(:,end),h,N,3,Ninterpolation);% Interpolate within each element to show internal oscillations
        subplot(2,2,1)
        plot(0:1/(length(Uinterpolated)):1,[Uinterpolated;Uinterpolated(1)],'b','Linewidth',3)
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(Uinterpolated)*L/nu)))
        xlabel('x/(2*\pi)'); ylabel('u(t)'); grid on
        xlim([0 1])

%        hold on
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,[0:(2/(N*Ninterpolation)):1],300)' ;
%        plot([0:(2/(N*Ninterpolation)):1]/2,sol_theory,'r')
      
%        plot([X;X(end)+h]/L,[u(1:3:3*N,end);u(1,end)],'r')
%        hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-Uinterpolated(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
              
        subplot(1,2,2)
        loglog(0:((length_vec+1)/2-1),spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics,'r','Linewidth',3, ...
               reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
        grid on; xlabel('k'); ylabel('E(k)')
        xlim([1 reference_spectrum(end,1)])

        drawnow
    end
    
    z=z+1;
    CFL=u(1:3:end,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        break;
    end
  end

%  relative_error

  spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  save(filename2,'spectralEnergyOut');
  
%  Uinterpolated = Interpolation(u(:,end),h,N,3,Ninterpolation);
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');

%  filename=strcat('Energie_',name,'.mat');
%  save(filename,'kinEnergy');
  
end

function y = RK4_FE_HermiteH5 (u,deltat,h,N,M,k,nu,F,constant_sub,ind)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
  
  Un=u;

%  w   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % weight for numerical integration
%    ksi = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % coordinate of point for numerical integration
%    for i=1:length(w)
%    d_shape_fct_vector(:,i) = get_deriv_shape_fct(ksi(i));    
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Cj = get_non_linear_hermite_H5(Un,ind,N,h) ;
%  Cj += constant_sub * h *h * get_subgrid_terms(h,Un2,ind,N,d_shape_fct_vector);
    
  k1 = - M \ (nu*k*Un + Cj);
  Un2 = Un + deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Cj = get_non_linear_hermite_H5(Un2,ind,N,h) ;
%  Cj += constant_sub * h *h * get_subgrid_terms(h,Un2,ind,N,d_shape_fct_vector);
  
  k2 = - M \ (nu*k*Un2 + Cj);
  Un3 = Un + deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  Cj = get_non_linear_hermite_H5(Un3,ind,N,h) ;
%  Cj += constant_sub * h *h * get_subgrid_terms(h,Un3,ind,N,d_shape_fct_vector);
  
  k3 = - M \ (nu*k*Un3 + Cj);
  Un4 = Un + deltat*k3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  Cj = get_non_linear_hermite_H5(Un4,ind,N,h) ;
%  Cj += constant_sub * h *h * get_subgrid_terms(h,Un4,ind,N,d_shape_fct_vector);
  
  k4 = - M \ (nu*k*Un4 + Cj);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 + F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cnonlinear = get_non_linear_hermite_H5(Un,ind,N,h)
% Compute the discretization of the nonlinear convection term u * (du/dx)
  Cnonlinear = zeros(3*N,1);
    
  u1 = Un(ind(:,1));  du1 = Un(ind(:,2));    ddu1 = Un(ind(:,3));
  u2 = Un(ind(:,4));  du2 = Un(ind(:,5));    ddu2 = Un(ind(:,6));
  u3 = Un(ind(:,7));  du3 = Un(ind(:,8));    ddu3 = Un(ind(:,9));
  
  Cnonlinear(1:3:3*N) = (-496*h^2*du1.^2 + 732*du1.*du2*h^2 - 732*du2.*du3*h^2 + 496*h^2*du3.^2 - ...
                         88*ddu1.*du1*h^3 - 72*ddu2.*du1*h^3 + 72*ddu1.*du2*h^3 + 176*ddu2.*du2*h^3 + ...
                         72*ddu3.*du2*h^3 - 72*ddu2.*du3*h^3 - 88*ddu3.*du3*h^3 - 4*h^4*ddu1.^2 - ...
                         7*ddu1.*ddu2*h^4 + 7*ddu2.*ddu3*h^4 + 4*h^4*ddu3.^2 - 3848*h*du1.*u1 + ...
                         2444*h*du2.*u1 - 328*h^2*ddu1.*u1 - 244*h^2*ddu2.*u1 - 8008*u1.^2 - 2444*h*du1.*u2 + ...
                         7696*h*du2.*u2 - 2444*h*du3.*u2 - 244*h^2*ddu1.*u2 + 244*h^2*ddu3.*u2 - 8008*u1.*u2 + ...
                         2444*h*du2.*u3 - 3848*h*du3.*u3 + 244*h^2*ddu2.*u3 + 328*h^2*ddu3.*u3 + 8008*u2.*u3 + ...
                         8008*u3.^2) /48048 ;
  Cnonlinear(2:3:3*N) = (2032*h^3*du1.^2 - 2032*du1.*du2*h^3 - 2032*du2.*du3*h^3 + 2032*h^3*du3.^2 + ...
                         345*ddu1.*du1*h^4 + 220*ddu2.*du1*h^4 - 184*ddu1.*du2*h^4 + 184*ddu3.*du2*h^4 - ...
                         220*ddu2.*du3*h^4 - 345*ddu3.*du3*h^4 + 15*ddu1.^2*h^5 + 20*ddu1.*ddu2*h^5 + ...
                         12*ddu2.^2*h^5 + 20*ddu2.*ddu3*h^5 + 15*ddu3.^2*h^5 + 16644*h^2*du1.*u1 - ...
                         7440*h^2*du2.*u1 + 1357*h^3*ddu1.*u1 + 805*h^3*ddu2.*u1 + 36660*h*u1.^2 + ...
                         5664*h^2*du1.*u2 - 5664*h^2*du3.*u2 + 502*h^3*ddu1.*u2 - 180*h^3*ddu2.*u2 + ...
                         502*h^3*ddu3.*u2 + 21060*h*u1.*u2 - 115440*h*u2.^2 + 7440*h^2*du2.*u3 - ...
                         16644*h^2*du3.*u3 + 805*h^3*ddu2.*u3 + 1357*h^3*ddu3.*u3 + 21060*h*u2.*u3 + ...
                         36660*h*u3.^2) /720720 ;
  Cnonlinear(3:3:3*N) = (-736*du1.^2*h^4 + 500*du1.*du2*h^4 - 500*du2.*du3*h^4 + 736*du3.^2*h^4 - 120*ddu1.*du1*h^5 - ...
                         60*ddu2.*du1*h^5 + 40*ddu1.*du2*h^5 - 48*ddu2.*du2*h^5 + 40*ddu3.*du2*h^5 - ...
                         60*ddu2.*du3*h^5 - 120*ddu3.*du3*h^5 - 5*ddu1.^2*h^6 - 5*ddu1.*ddu2*h^6 + 5*ddu2.*ddu3*h^6 + ...
                         5*ddu3.^2*h^6 - 6328*h^3*du1.*u1 + 2060*h^3*du2.*u1 - 496*h^4*ddu1.*u1 - 240*h^4*ddu2.*u1 - ...
                         14640*h^2*u1.^2 - 1108*h^3*du1.*u2 - 9840*h^3*du2.*u2 - 1108*h^3*du3.*u2 - 76*h^4*ddu1.*u2 + ...
                         76*h^4*ddu3.*u2 - 5040*h^2*u1.*u2 + 2060*h^3*du2.*u3 - 6328*h^3*du3.*u3 + 240*h^4*ddu2.*u3 + ...
                         496*h^4*ddu3.*u3 +  5040*h^2*u2.*u3 + 14640*h^2*u3.^2) /2882880 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%