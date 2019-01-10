function FD_compact_spectral (N,nu,constant_sub,L,time,nbrpointtemp,name,file_spectrum)
% Solve the 1D forced Burgers equation with compact finite difference schemes
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp('*********************************************************************')  
  disp('Finite difference - Compact schemes with spectral-like resolution')
  disp('*********************************************************************')  

  h = L/(N);% Length of the elements
  X = linspace(0,L,N)';
  nbrpointtime = nbrpointtemp;
  deltat           = time / nbrpointtime;

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1) = 2*rand(N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(:,1) = sin(X * 2*pi/L);

  timeBeforeStatistics = 10;

% Multiplicative factors for the first and second spatial derivatives
  factors_deriv = zeros(5,2) ;
% Optimal - see the doc
  factors_deriv(1:5,1) = [0.5771439  0.0896406  1.3025166  0.99355  0.03750245]'; % alpha beta a b c
  factors_deriv(1:5,2) = [0.50209266 0.05569169 0.21564935 1.723322 0.17659730]'; % alpha beta a b c
% Order 10 - see the doc
%  factors_deriv(1:5,1) = [0.5  1/20  17/12  101/150  0.01]';
%  factors_deriv(1:5,2) = [334/899 43/1798 1065/1798 1038/899 79/1798]';
% Order 8 - see the doc
%  factors_deriv(1:5,1) = [4/9 1/36 40/27 25/54 0]';
%  alpha = 344/1179;
%  factors_deriv(1:5,2) = [alpha (38*alpha-9)/214 (696-1191*alpha)/428 (2454*alpha-294)/535 0]';
% Order 6 - see the doc
%  alpha = 1/3;  beta=0;
%  factors_deriv(1:5,1) = [alpha beta (9+alpha-20*beta)/6 (-9+32*alpha+62*beta)/15 (1-3*alpha+12*beta)/10]';
%  alpha = 2/11; beta=0;
%  factors_deriv(1:5,2) = [alpha beta (6-9*alpha-12*beta)*0.25 (-3+24*alpha-6*beta)*0.2 (2-11*alpha+124*beta)/20]';
% Order 4 central scheme - see the doc
%  alpha = 0;  beta=0;
%  factors_deriv(1:5,1) = [0 0 2*(2+alpha)/3 (4*alpha-1)/3  0]';
%  alpha = 0; beta=0;
%  factors_deriv(1:5,2) = [0 0 4*(1-alpha)/3 (10*alpha-1)/3  0]';
% Order 2 central scheme - see the doc
%  factors_deriv(1:5,1) = [0 0 1 0  0]';
%  factors_deriv(1:5,2) = [0 0 1 0  0]';

% Scale the factors
  factors_deriv(3,1) = factors_deriv(3,1)/(2*h); factors_deriv(4,1) = factors_deriv(4,1)/(4*h);   factors_deriv(5,1) = factors_deriv(5,1)/(6*h); 
  factors_deriv(3,2) = factors_deriv(3,2)/(h*h); factors_deriv(4,2) = factors_deriv(4,2)/(4*h*h); factors_deriv(5,2) = factors_deriv(5,2)/(9*h*h); 
  factors_deriv(3:5,2) = factors_deriv(3:5,2) * nu ;
  
% Create the matrix for the first spatial derivative
  diag0 = ones(N,1);
  diag1 = factors_deriv(1,1) * ones(N,1); % alpha
  diag2 = factors_deriv(2,1) * ones(N,1); % beta
  C0 = diag(sparse(diag0));  
  Cp1 = diag(sparse(diag1),1);  Cm1 = diag(sparse(diag1),-1);
  Cp2 = diag(sparse(diag2),2);  Cm2 = diag(sparse(diag2),-2);
  mat_deriv1 = C0 + Cp1(1:N,1:N) + Cm1(2:N+1,2:N+1) + Cp2(1:N,1:N) + Cm2(3:N+2,3:N+2);
  mat_deriv1(N-1,1) = diag2(N-1) ;   mat_deriv1(2,N)     = diag2(2) ;
  mat_deriv1(N,1)   = diag1(N);      mat_deriv1(N,2)     = diag2(N);
  mat_deriv1(1,N)   = diag1(1);      mat_deriv1(1,N-1)   = diag2(1); 
  
% Create the matrix for the second spatial derivative
  diag0 = ones(N,1);
  diag1 = factors_deriv(1,2) * ones(N,1); % alpha
  diag2 = factors_deriv(2,2) * ones(N,1); % beta
  C0 = diag(sparse(diag0));  
  Cp1 = diag(sparse(diag1),1);  Cm1 = diag(sparse(diag1),-1);
  Cp2 = diag(sparse(diag2),2);  Cm2 = diag(sparse(diag2),-2);
  mat_deriv2 = C0 + Cp1(1:N,1:N) + Cm1(2:N+1,2:N+1) + Cp2(1:N,1:N) + Cm2(3:N+2,3:N+2);
  mat_deriv2(N-1,1) = diag2(N-1) ;   mat_deriv2(2,N)     = diag2(2) ;
  mat_deriv2(N,1)   = diag1(N);      mat_deriv2(N,2)     = diag2(N);
  mat_deriv2(1,N)   = diag1(1);      mat_deriv2(1,N-1)   = diag2(1);
      
  kinEnergy = zeros(nbrpointtemp+1,1); kinEnergy(1) = u(:,1)' * u(:,1) * h * 0.5 ;
% energy_conv is the numerical energy produced by the spatial discretization of the convective term,
% it must be equal to zero for an energy conservative scheme
  energy_conv= zeros(nbrpointtemp+1,1);
    
%  filename=strcat(name,num2str(1),'.mat'); uu=u(:,1); save(filename,'uu');

  z=2; j=2; ind_error=1;
  nbrPointsStatistics=0;  kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);

  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand();    phi3=2*pi*rand();
    KForce=2*pi/L;
    F = sqrt(deltat)*(cos(2*KForce*X+phi2)+cos(3*KForce*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;
    
%******** Call Runge-Kutta and compute kinematic energy ********
    [u(:,z), energy_conv(i)] = RK4_FD_compact_spectral (u(:,z-1),deltat,N,mat_deriv1,mat_deriv2,F,h,constant_sub,factors_deriv);
    kinEnergy(i) = h*0.5*u(:,z)'*u(:,z);
        
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);
      fft_u = fft(u(:,end));
      newSpectralEnergy = fft_u.*conj(fft_u)/N/N;
      spectralEnergy = spectralEnergy + 2*pi*newSpectralEnergy;
      nbrPointsStatistics = nbrPointsStatistics + 1;
    end
    
% Save the results, free some memory and show the results
    if ( mod(z/nbrpointtime,0.1) == 0)
        uu=u(:,end);
%        filename=[name,num2str(j),'.mat']; save(filename,'uu'); j=j+1;
        clear u;
        z=1;
        u(:,1)=uu;
        
        reynolds_number        = mean(uu)*L/nu;
        kolmogorov_length      = L*reynolds_number^(-0.75);
        kolmogorov_wave_length = 2*pi/kolmogorov_length ;
                
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ));
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'-ascii','spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ));
        end
        
        subplot(2,2,1);
        plot(X/L, uu)
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(uu)*L/nu)))
        
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-uu(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
        
        subplot(2,2,3)
        plot((0:(i-1))*deltat,kinEnergy(1:i))
        grid on; xlabel('Time'); ylabel('E(t)')
         
        subplot(2,2,2)
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics, reference_spectrum(:,1),reference_spectrum(:,2))
        grid on; xlabel('k'); ylabel('E(k)')
        
        subplot(2,2,4)
        plot((0:(i-1))*deltat,energy_conv(1:i),'b')
        grid on; xlabel('Time'); ylabel('E_{prod}(t)')
            
        drawnow;
    end
    
    z=z+1;
  
    CFL=u(:,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(strcat('Divergence of ',name));
        break;
    end
    end
      
%  relative_error
  
   spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
   filename2=strcat('Spectral_energy_',name,'.mat');
   save(filename2,'-ascii','spectralEnergyOut');
  
  %filename=strcat('Energy_',name,'.mat');
  %save(filename,'kinEnergy');
  
end
