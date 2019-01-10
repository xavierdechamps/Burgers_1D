function FE_HermiteP3(N,nu,constant_sub,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
% Solve the 1D forced Burgers equation with cubic Hermite elements 
% The unknowns are the velocity and the first spatial derivative of the velocity, thus 2 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element cubic Hermite P3")
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
    u(:,z) = RK4_FE_HermiteP3 (u(:,z-1),deltat,h,N,M,K,nu,F,constant_sub);
    
    kinEnergy(i)=get_kinematic_energy(h,DG,u(:,end),N,2);
    
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
        subplot(1,2,1)
        plot(0:1/(length(Uinterpolated)):1,[Uinterpolated;Uinterpolated(1)],'b','Linewidth',3)
        title(strcat('Time= ',num2str(i*deltat),', Re= ',num2str(mean(Uinterpolated)*L/nu)))
        xlabel('x/(2*\pi)'); ylabel('u(t)'); grid on
          
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
        
        subplot(1,2,2) ;
        loglog(0:N-1,spectralEnergy(1:N)/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
        grid on; xlabel('k'); ylabel('E(k)');

        drawnow ;
    end
    
    z=z+1;
    CFL=u(1:2:end,end).*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(strcat('Divergence of ',name));
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
