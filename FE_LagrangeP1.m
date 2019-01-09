function FE_LagrangeP1(N,nu,constant_sub,L,time,nbrpointtemp,name,file_spectrum)
% Solve the 1D forced Burgers equation with cubic Hermite elements 
% The unknowns are the velocity and the first spatial derivative of the velocity, thus 2 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element linear Lagrange P1")
  disp("************************************************************")

  h=L/N;% Length of the elements
  DG=2;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;
  
% ******** Mass Mij and rigidity Kij matrices *************
  Mij=[2 1;1 2].*h/6;
  Kij=[1 -1;-1 1]./h;
  
% ************* Assemble the matrices (with periodic condition) ***********
  M=sparse(N,N);
  k=sparse(N,N);
  for i=1:N-1
    M(i:i+1,i:i+1)=M(i:i+1,i:i+1)+Mij;
    k(i:i+1,i:i+1)=k(i:i+1,i:i+1)+Kij;
  end
  M(1,1)=M(1,1)+Mij(2,2); M(N,N)=M(N,N)+Mij(1,1); M(1,N)=M(1,N)+Mij(2,1); M(N,1)=M(N,1)+Mij(1,2);
  k(1,1)=k(1,1)+Kij(2,2); k(N,N)=k(N,N)+Kij(1,1); k(1,N)=k(1,N)+Kij(2,1); k(N,1)=k(N,1)+Kij(1,2);

%%%%%%%%%%%%%%%%%%   WARNING %%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%   LUMP MASS MATRIX USED TO REDUCE THE TIME INCREMENT
%%%%%%%%%%%%%%%%%%   THE DISCRETIZATION IS THEN IDENTICAL TO THE CONSERVATIVE FINITE DIFFERENCE HC2
%  M = eye(N) * h ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1)=2*rand(N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(:,1)=sin(X * 2*pi/L);
%  u(:,1)=2*exp(-(X-3).^2);
%  u(:,1) = sin(pi*X) 

  timeBeforeStatistics = 10;

  j=2;   z=2;  ind_error = 1;
  filename=[name,num2str(1),'.mat'];
  uu=u(:,1);
  %save(filename,'uu');
  
  kinEnergy     = zeros(1:nbrpointtime+1,1);
  kinEnergy(1) = get_kinematic_energy(h,DG,u(:,1),N,1);
  nbrPointsStatistics=0;  kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);

  %[maxU, maxInd] = max(u(:,1)); [minU, minInd] = min(u(:,1));
  %distance_sinus = zeros(1:nbrpointtime+1,1);
  %distance_sinus(1,1:4) = [0 maxU minU X(minInd)-X(maxInd)];

  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2   = 2*pi*rand();    phi3   = 2*pi*rand();
    KForce = 2*pi/L;
    F      = sqrt(deltat) * ( cos(2*KForce*X+phi2) + cos(3*KForce*X+phi3) );
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;
        
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z)     = RK4_FE_Lagrangep1(u(:,z-1),deltat,N,M,k,nu,F,constant_sub);
    
    kinEnergy(i) = get_kinematic_energy(h,DG,u(:,z),N,1);
    
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);
      fft_u = fft(u(:,end));
      nlleEnergieSpect = fft_u.*conj(fft_u)/N/N;
      spectralEnergy = spectralEnergy + 2*pi*nlleEnergieSpect;
      nbrPointsStatistics = nbrPointsStatistics + 1;
    end
    
  %    [maxU, maxInd] = max(u(:,z)); [minU, minInd] = min(u(:,z));
  %    distance_sinus(i,1:4) = [(i-1)*deltat maxU minU X(minInd)-X(maxInd)];
    
% Save the results, free some memory and show the results
    if ( mod( (i-1)*deltat , 0.1) == 0)
        uu=u(:,end);
%        filename=[name,num2str(j),'.mat']; save(filename,'uu'); j=j+1;
        clear u;
        z=1;
        u(:,1)=uu;
        
%        reynolds_number        = mean(uu)*L/nu;
%        kolmogorov_length      = L*reynolds_number^(-0.75);
%        kolmogorov_wave_length = 2*pi/kolmogorov_length
        
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=['Spectral_energy_',name,'.mat']; energieSpectraleOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'energieSpectraleOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
                
%        disp(strcat("Max at x=",num2str(X(maxInd)), " and min at x=",num2str(X(minInd))," and distance=",num2str(distance_sinus(i,4))));
        
        subplot(2,2,1)
        plot([X; X(end)+h]/L,[uu; uu(1)])
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(uu)*L/nu)))
                
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-uu(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
      
        subplot(2,2,3);
        plot((1:i)*deltat,kinEnergy(1:i));
        grid on; xlabel('Time'); ylabel('E(t)');
                        
        subplot(2,2,2);
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics, reference_spectrum(:,1),reference_spectrum(:,2));
        grid on; xlabel('k'); ylabel('E(k)');
        
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

  %relative_error

  energieSpectraleOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Energie_Spectrale_',name,'.mat');
  save(filename2,'energieSpectraleOut');

%  tosave = [X u(:,end)];
%  save(strcat(name,'.mat'),'tosave');
        
%  save(strcat(name,'_distance_sinus.mat'),'distance_sinus');

%  filename=strcat('Energie_',name,'.mat');
%  save(filename,'kinEnergy');
end
