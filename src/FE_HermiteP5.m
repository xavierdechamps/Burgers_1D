function FE_HermiteP5(N,nu,constant_sub,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
% Solve the 1D forced Burgers equation with 5th order Hermite elements 
% The unknowns are the velocity, the first and second spatial derivatives of the velocity,
% thus 3 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element 5th order Hermite P5")
  disp("************************************************************")

  length_vec = 3*N;
  h  = L/N; % Length of the elements
  DG = 4;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;

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
    u(:,z) = RK4_FE_HermiteP5 (u(:,z-1),deltat,h,N,M,K,nu,F,constant_sub);
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
        loglog(0:N-1,spectralEnergy(1:N)/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
        grid on; xlabel('k'); ylabel('E(k)')

        drawnow
    end
    
    z=z+1;
    CFL=u(1:3:end,end).*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(strcat('Divergence of ',name));
        break;
    end
  end

%  relative_error

  spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics;
  filename2=strcat('Energie_Spectrale_',name,'.mat');
  save(filename2,'spectralEnergyOut');
  
%  Uinterpolated = Interpolation(u(:,end),h,N,3,Ninterpolation);
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');

%  filename=strcat('Energie_',name,'.mat');
%  save(filename,'kinEnergy');
  
end
