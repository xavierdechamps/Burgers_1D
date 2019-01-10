function FE_LagrangeP3(N,nu,constant_sub,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
% Solve the 1D forced Burgers equation with cubic Hermite elements 
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element cubic Lagrange P3")
  disp("************************************************************")
  
  length_vec = 3*N;
  h = L/N;% Length of the elements
  DG=4;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,3*N)';

  nbrpointtime = nbrpointtemp;
  deltat       = time / nbrpointtime;

% ******** Mass Mij and rigidity Kij matrices *************
  Mij = h * [8/105    33/560  -3/140    19/1680  ;
            33/560   27/70   -27/560  -3/140    ;
           -3/140   -27/560   27/70    33/560   ;
            19/1680 -3/140    33/560   8/105     ];

% LUMPED MATRIX
%  Mij = h * [210/1680     0         0          0     ;
%              0      630/1680     0          0     ;
%              0         0      630/1680      0     ;
%              0         0         0       210/1680 ];
              
              
  Kij = [ 37/10   -189/40   27/20   -13/40  ;
        -189/40   54/5    -297/40   27/20  ;
         27/20   -297/40   54/5    -189/40 ;
        -13/40    27/20   -189/40   37/10        ]/h;
       
% ************* Assemble the matrices (with periodic condition) ***********
  M = sparse(length_vec,length_vec);
  K = sparse(length_vec,length_vec);
  for i=1:N
    indmin = 3*i-2;
    indmax = 3*i+1;
    ind    = 4;
    if (i==(N))
      indmax = 3*i; 
      ind = 3;
    end
    M(indmin:indmax,indmin:indmax) = M(indmin:indmax,indmin:indmax) + Mij(1:ind,1:ind);
    K(indmin:indmax,indmin:indmax) = K(indmin:indmax,indmin:indmax) + Kij(1:ind,1:ind);
  end
  % last element for periodicity
  M(1,end-2)=M(1,end-2)+Mij(4,1);  M(1,end-1)=M(1,end-1)+Mij(4,2);  M(1,end)=M(1,end)+Mij(4,3);  M(1,1)=M(1,1)+Mij(4,4);
  M(end,1)  =M(end,1)  +Mij(3,4);  M(end-1,1)=M(end-1,1)+Mij(2,4);  M(end-2,1)=M(end-2,1)+Mij(1,4);

  K(1,end-2)=K(1,end-2)+Kij(4,1);  K(1,end-1)=K(1,end-1)+Kij(4,2);  K(1,end)=K(1,end)+Kij(4,3);  K(1,1)=K(1,1)+Kij(4,4);
  K(end,1)  =K(end,1)  +Kij(3,4);  K(end-1,1)=K(end-1,1)+Kij(2,4);  K(end-2,1)=K(end-2,1)+Kij(1,4);

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u = 2*rand(length_vec,1) - 1 ;
% Sinus solution for non-forced Burgers equation
%  u(:,1)=sin(X(:) * 2*pi/L);
%  u(:,1)=2*exp(-(X(:)-3).^2);

  timeBeforeStatistics = 10;

% kinematic energy not computed yet for the cubic Lagrange element
%  kinEnergy    = zeros(1:nbrpointtime+1,1);  kinEnergy(1) = get_kinematic_energy(h,DG,u(:,1),3*N,1);
  j=2;   z=2;  ind_error=1;
  filename=[name,num2str(1),'.mat'];
  uu=u(:,1);
%  save(filename,'uu');

  nbrPointsStatistics=0;

  spectralEnergy=zeros(Ninterpolation*N,1);
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
%      F = 0;
        
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z)     = RK4_FE_Lagrangep3 (u(:,z-1),deltat,h,N,M,K,nu,F,constant_sub);
    
% kinematic energy not computed yet for the cubic Lagrange element
%      kinEnergy(i) = get_kinematic_energy(h,DG,u(:,z),3*N,1);
    
    if (i*deltat>=timeBeforeStatistics)
      Uinterpolated        = Interpolation(u(:,end),h,N,2,Ninterpolation);% Interpolate within each element to show internal oscillations
      length_Uinterpolated = length(Uinterpolated);
      fft_u                = fft(Uinterpolated);
      newSpectralEnergy    = fft_u.*conj(fft_u)/length_Uinterpolated/length_Uinterpolated;
      spectralEnergy       = spectralEnergy + 2*pi*newSpectralEnergy;
      nbrPointsStatistics  = nbrPointsStatistics + 1;
    end
    
%    [maxU, maxInd] = max(u(:,z)); [minU, minInd] = min(u(:,z));
%    distance_sinus(i,1:4) = [(i-1)*deltat maxU minU X(minInd)-X(maxInd)];
    
% Save the results, free some memory and show the results
    if ( mod(z/nbrpointtime,0.1) == 0)
        uu=u(:,end);
%        filename=[name,num2str(j),'.mat']; save(filename,'uu'); j=j+1;
        clear u;
        z=1;
        u(:,1)=uu;
        
%        reynolds_number        = mean(uu)*L/nu;
%        kolmogorov_length      = L*reynolds_number^(-0.75);
%        kolmogorov_wave_length = 2*pi/kolmogorov_length
        
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ));
%           filename2=['Spectral_energy_',name,'.mat']; spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ));
        end
                
          subplot(2,2,1);
          Uinterpolated = Interpolation(u(:,end),h,N,2,Ninterpolation);% Interpolate within each element to show internal oscillations
          plot(0:1/(length(Uinterpolated)):1,[Uinterpolated;Uinterpolated(1)],'b','Linewidth',3);
          grid on; xlabel('x/(2*\pi)'); ylabel('u(t)');
          title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(uu)*L/nu)));
%          hold on
%          plot([X;X(end)+h/3]/L,[u(:,end); u(1,end)],'b*')
%          hold off
       
%          sol_theory = get_analytical_solution((i-1)*deltat,nu,[0:(2/(N*Ninterpolation)):1],100)' ;
%          plot([0:(2/(N*Ninterpolation)):1]/2,sol_theory,'r')
%          hold off
%          relative_error(ind_error,1) = (i-1)*deltat;
%          relative_error(ind_error,2) = sqrt( sum( (sol_theory-Uinterpolated(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%          disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%          ind_error = ind_error + 1;
        
%          subplot(2,2,3)
%          plot((1:i)*deltat,kinEnergy(1:i))
%          grid on; xlabel('Time'); ylabel('E(t)')
         
          subplot(1,2,2);
          loglog(0:((length_vec+1)/2-1),spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
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

%  relative_error

  spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics;
  filename2=['Spectral_energy_',name,'.mat'];
  save(filename2,'spectralEnergyOut');

%  Uinterpolated = Interpolation(u(:,end),h,N,2,Ninterpolation);%Raffinage de l'affichage
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');
        
%  filename=strcat('Energie_',name,'.mat');
%  save(filename,'kinEnergy');
  
end
