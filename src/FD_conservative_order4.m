function FD_conservative_order4 (N,nu,constant_sub,L,time,nbrpointtemp,name,file_spectrum)
% Solve the 1D forced Burgers equation with the energy conservative Hs scheme of 
% order 4 for the convective term
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("*********************************************************************")  
  disp("Finite difference - 4th order conservative scheme for convective term")
  disp("*********************************************************************")  

  h = L/N;% Length of the elements
  X = linspace(0,L,N)';
  nbrpointtime = nbrpointtemp;
  deltat       = time / nbrpointtime;

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1) = 2*rand(N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(:,1)=sin(X(:)* 2*pi/L);

  timeBeforeStatistics = 10;

% ******** Rigidity Kij matrix *************
  K=sparse(N,N);
  for i=1:N 
% Second order
%    switch i
%      case 1
%        a=[N 1 2];
%      case N
%        a=[N-1 N 1];
%      otherwise
%        a=[i-1 i i+1];
%    end
%    K(i,a) = [1 -2 1]*nu/(h*h);
  
% Fourth order
    switch i
      case 1
        a=[N-1 N 1 2 3]; % periodic condition on the first node
      case 2
        a=[N 1 2 3 4];% periodic condition on the second node
      case N-1
        a=[N-3 N-2 N-1 N 1];
      case N
        a=[N-2 N-1 N 1 2];% periodic condition on the last node
      otherwise
        a=[i-2 i-1 i i+1 i+2]; % indices of 5 successive nodes
    end
    K(i,a) = [-1 16 -30 16 -1]*nu/(12*h*h);

% sixth order
%    switch i
%      case 1
%        a=[N-2 N-1 N 1 2 3 4];
%      case 2
%        a=[N-1 N 1 2 3 4 5];
%      case 3
%        a=[N 1 2 3 4 5 6];
%      case N-2
%        a=[N-5 N-4 N-3 N-2 N-1 N 1];
%      case N-1
%        a=[N-4 N-3 N-2 N-1 N 1 2];
%      case N
%        a=[N-3 N-2 N-1 N 1 2 3];
%      otherwise
%        a=[i-3 i-2 i-1 i i+1 i+2 i+3];
%    end
%    K(i,a) = [1/90   -3/20   1.5   -49/18   1.5   -3/20   1/90]*nu/(h*h);
  end
  
  kinEnergy = zeros(nbrpointtemp+1,1); kinEnergy(1) = h*0.5*u(:,1)'*u(:,1);
  
%  filename=[name,num2str(1),'.mat']; uu=u(:,1); save(filename,'uu');

  z=2; j=2; ind_error=1;
  nbrPointsStatistics=0; kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);

  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand();    phi3=2*pi*rand();
    KForce=2*pi/L;
    F = sqrt(deltat)*(cos(2*KForce.*X+phi2)+cos(3*KForce.*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F=0;
    
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z) = RK4_FD_conservative_order4 (u(:,z-1),deltat,N,K,F,h,constant_sub);

    kinEnergy(i) = h*0.5*u(:,z)'*u(:,z);
        
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);
      fft_u = fft(u(:,end));
      newSpectralEnergy = fft_u.*conj(fft_u)/(N)/(N);
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
        
%        reynolds_number        = mean(uu)*L/nu;
%        kolmogorov_length      = L*reynolds_number^(-0.75);
%        kolmogorov_wave_length = 2*pi/kolmogorov_length ;
        
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
                
        subplot(2,2,1)
        plot([X;L]/L, [uu; uu(1)],'Linewidth',3)
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        title(strcat('Time= ',num2str(i*deltat),', Re= ',num2str(mean(uu)*L/nu)))
        
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-uu(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
        
         subplot(2,2,3)
         plot((0:(i-1))*deltat,kinEnergy(1:i),'Linewidth',3)
         grid on; xlabel('Time'); ylabel('E(t)')
         
         subplot(1,2,2)
         loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
         grid on; xlabel('k'); ylabel('E(k)')
        
         drawnow
    end
    
    z=z+1;
  
    CFL=u(:,end).*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(strcat('Divergence of ',name));
        break;
    end
  end
  
%  relative_error
  
  spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  save(filename2,'spectralEnergyOut');

  %figure()
  %loglog(spectralEnergyOut)
  %grid on

  %filename=strcat('Energie_',name,'.mat');
  %save(filename,'kinEnergy');
  
end
