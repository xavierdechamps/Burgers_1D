function FE_LagrangeP1(N,nu,constant_sub,filter,Alpha_Pade,L,time,nbrpointtemp,name,file_spectrum)
% Solve the 1D forced Burgers equation with cubic Hermite elements 
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element linear Lagrange P1")
  switch filter
     case 0
       disp("   Constant value Smagorinsky model")
     case 1
       disp("   Dynamic Smagorinsky model - 3 points stencil for the low-pass filter")
     case 2
       disp("   Dynamic Smagorinsky model - 5 points stencil for the low-pass filter")
     case 3
       disp("   Dynamic Smagorinsky model - 7 points stencil for the low-pass filter")
     case 4
       disp("   Dynamic Smagorinsky model - 9 points stencil for the low-pass filter")
     case 5
       disp("   Dynamic Smagorinsky model - Pade low-pass filter")
     otherwise 
       disp("   Direct numerical simulation")
  end
  disp("************************************************************")

  h=L/N;% Length of the elements
  DG=2;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;
  
% ************* Assemble the matrix (with periodic condition) ***********
  Mij =     [2 1;1 2]*h/6;
  M   = sparse(N,N);
  for i=1:N-1
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + Mij;
  end
  M(1,1) += Mij(2,2); M(N,N) += Mij(1,1); M(1,N) += Mij(2,1); M(N,1) += Mij(1,2);

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
%  save(filename,'uu');
  
  kinEnergy    = zeros(1:nbrpointtime+1,1);
  kinEnergy(1) = get_kinematic_energy(h,DG,u(:,1),N,1);
  nbrPointsStatistics=0;  kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);

% Store the indices of the neighbour nodes i-4 i-3 i-2 i-1 i i+1 i+2 i+3 i+4 used in the derivatives
  ind = zeros(N,9);
  ind(:,5) = 1:N; % i
  ind(:,4) = circshift(ind(:,5),1,1) ; % i-1
  ind(:,3) = circshift(ind(:,5),2,1) ; % i-2
  ind(:,2) = circshift(ind(:,5),3,1) ; % i-3
  ind(:,1) = circshift(ind(:,5),4,1) ; % i-4
  ind(:,6) = circshift(ind(:,5),-1,1); % i+1
  ind(:,7) = circshift(ind(:,5),-2,1); % i+2
  ind(:,8) = circshift(ind(:,5),-3,1); % i+3
  ind(:,9) = circshift(ind(:,5),-4,1); % i+4
  dynamic_smag_constant = zeros(nbrpointtemp,1);
  mat_alpha = zeros(N,N) ;
  for i=1:N
     mat_alpha(i, ind(i,4:6)) = [Alpha_Pade , 1 , Alpha_Pade] ;
  end
  mat_alpha = sparse(mat_alpha);
  
  %[maxU, maxInd] = max(u(:,1)); [minU, minInd] = min(u(:,1));
  %distance_sinus = zeros(1:nbrpointtime+1,1);
  %distance_sinus(1,1:4) = [0 maxU minU X(minInd)-X(maxInd)];

  diverged = false;
  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2   = 2*pi*rand();    phi3   = 2*pi*rand();
    KForce = 2*pi/L;
    F      = sqrt(deltat) * ( cos(2*KForce*X+phi2) + cos(3*KForce*X+phi3) );
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;
        
%******** Call Runge-Kutta and compute kinematic energy ********
    [u(:,z),dynamic_smag_constant(i-1)] = ...
              RK4_FE_Lagrangep1(u(:,z-1),deltat,N,M,nu,h,F,constant_sub,ind,filter,Alpha_Pade ,mat_alpha);
    
    kinEnergy(i) = get_kinematic_energy(h,DG,u(:,z),N,1);
    
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);
      fft_u = fft(u(:,end));
      newSpectralEnergy = fft_u.*conj(fft_u)/N/N;
      spectralEnergy = spectralEnergy + 2*pi*newSpectralEnergy;
      nbrPointsStatistics = nbrPointsStatistics + 1;
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
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=['Spectral_energy_',name,'.mat']; spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
                
%        disp(strcat("Max at x=",num2str(X(maxInd)), " and min at x=",num2str(X(minInd))," and distance=",num2str(distance_sinus(i,4))));
        
        subplot(2,2,1)
        plot([X; X(end)+h]/L,[uu; uu(1)],'b','Linewidth',3)
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        xlim([0 1])
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(uu)*L/nu)))
                
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-uu(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
      
        subplot(2,2,3);
        plot((1:i)*deltat,kinEnergy(1:i),'b','Linewidth',3);
        grid on; xlabel('Time'); ylabel('E(t)');
        xlim([0 time])
                        
        subplot(2,2,2);
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
        grid on; xlabel('k'); ylabel('E(k)');
        xlim([1 reference_spectrum(end,1)])
        
        subplot(2,2,4)
        mean_Smagorinsky = mean(dynamic_smag_constant(1:i-1));
        standard_deviation = std(dynamic_smag_constant(1:i-1));
        standard_deviationp = mean_Smagorinsky + standard_deviation;
        standard_deviationm = mean_Smagorinsky - standard_deviation;
        plot((1:(i-1))*deltat,dynamic_smag_constant(1:i-1),'b','Linewidth',3) ; hold on;
        plot([1 (i-1)]*deltat,[mean_Smagorinsky mean_Smagorinsky],      'r-', 'Linewidth',3);
        plot([1 (i-1)]*deltat,[standard_deviationp standard_deviationp],'r--','Linewidth',3);
        plot([1 (i-1)]*deltat,[standard_deviationm standard_deviationm],'r--','Linewidth',3);
        hold off;
        grid on; xlabel('Time'); ylabel('Smagorinsky C_s(t)') ;
        xlim([0 time])
        
        drawnow;
    end
    z=z+1;
  
    CFL=u(:,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        diverged = true;
        break;
    end
  end

  mean_Smagorinsky = mean(dynamic_smag_constant)
  standard_deviation = std(dynamic_smag_constant)
  %relative_error

  spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  if (~diverged)
    save(filename2,'-ascii','spectralEnergyOut');
  end

%  tosave = [X u(:,end)];
%  save(strcat(name,'.mat'),'tosave');
        
%  save(strcat(name,'_distance_sinus.mat'),'distance_sinus');

%  filename=strcat('Energy_',name,'.mat');
%  save(filename,'kinEnergy');
end

function [y,smag_sub] = RK4_FE_Lagrangep1 (u,deltat,N,M,nu,h,F,constant_sub,ind,filter,alpha,mat_alpha)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with linear Lagrange elements
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
  Un = u;
  nu_over_h = nu / h ;
  
%%%%%% Get the Smagorinsky constant in case of dynamic model
  if (filter>0)
     kappa = 2; % filter ratio
     smag_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,alpha,mat_alpha);
  elseif (filter==0)
     smag_sub = constant_sub ;
  else 
     smag_sub = 0. ;
  end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Kj = get_viscous_term  (Un,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un,ind,smag_sub) ;
  k1 = M \ ( Kj + Cj );

  Un2 = Un + deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Kj = get_viscous_term  (Un2,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un2,ind,smag_sub) ;
  k2 = M \ ( Kj + Cj );
  
  Un3 = Un + deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  Kj = get_viscous_term  (Un3,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un3,ind,smag_sub) ;
  k3 = M \ ( Kj + Cj );

  Un4 = Un + deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  Kj = get_viscous_term  (Un4,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un4,ind,smag_sub) ;
  k4 = M \ ( Kj + Cj );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 + F;
end

function vecC = get_viscous_term(Un,ind,nu_over_h)
% nu d2(u)/d(x)2 
% ind has size (N,3)
   vecC = nu_over_h * ( Un(ind(:,1)) - 2*Un(ind(:,2)) + Un(ind(:,3)) ) ;
endfunction

function vecC = get_nonlinear_term(Un,ind,constant_sub)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
% Convective term - u du/dx
   vecC = ( Un(ind(:,4)) .* ( Un(ind(:,5)) + Un(ind(:,4)) ) - Un(ind(:,6)) .* ( Un(ind(:,5)) + Un(ind(:,6)) ) ) / 6 ;
   
% Subgrid term
   if ( constant_sub>0 )
      vecC += constant_sub^2 * ( abs(Un(ind(:,6))-Un(ind(:,5))).*(Un(ind(:,6))-Un(ind(:,5))) - ...
                                 abs(Un(ind(:,5))-Un(ind(:,4))).*(Un(ind(:,5))-Un(ind(:,4))) ) ;
   endif
endfunction

function dudx = get_first_derivative(Un,ind,h)
% ind has size (N,3)
   dudx = ( Un(ind(:,3)) - Un(ind(:,1)) ) * 0.5 / h ;
endfunction

function smooth = apply_filter(Un,ind,type,alpha,mat_alpha)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
   N = length(Un);
   switch type
      case 1
% Low-pass filter binomial over 3 points B2
         smooth = 0.25 * ( Un(ind(:,4)) + 2*Un(ind(:,5)) + Un(ind(:,6)) ) ;
      case 2
% Low-pass filter binomial over 5 points B(2,1)
         smooth = ( -Un(ind(:,3)) + 4*Un(ind(:,4)) + 10*Un(ind(:,5)) + 4*Un(ind(:,6)) - Un(ind(:,7)) )/16;
      case 3
% Low-pass filter binomial over 7 points B(3,1)
         smooth = ( Un(ind(:,2)) - 6*Un(ind(:,3)) + 15*Un(ind(:,4)) + 44*Un(ind(:,5)) + ...
                            15*Un(ind(:,6)) - 6*Un(ind(:,7)) + Un(ind(:,8)) )/64;
      case 4
% Low-pass filter binomial over 9 points B(4,1)
         smooth = ( -Un(ind(:,1)) +8*Un(ind(:,2)) - 28*Un(ind(:,3)) + 56*Un(ind(:,4)) + 186*Un(ind(:,5)) + ...
                     56*Un(ind(:,6)) - 28*Un(ind(:,7)) + 8*Un(ind(:,8)) - Un(ind(:,9)) )/256;
     case 5
% Pade filter
         a0 = (11 + 10*alpha)/32;
         a1 = (15 + 34*alpha)/64;
         a2 = (-3 + 6*alpha)/32;
         a3 = ( 1 - 2*alpha)/64;
         RHS = 2*a0*Un(ind(:,5))             + a1*(Un(ind(:,4))+Un(ind(:,6))) + ...
              a2*(Un(ind(:,3))+Un(ind(:,7))) + a3*(Un(ind(:,2))+Un(ind(:,8))) ;
         smooth = mat_alpha \ RHS ;
      otherwise
          disp("Unknown type of filter");
          smooth = Un ;
   end
endfunction

function dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,alpha,mat_alpha)
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
   u_filter = apply_filter(Un ,ind,filter,alpha,mat_alpha) ;
   L        = apply_filter(Un.*Un ,ind,filter,alpha,mat_alpha) - u_filter.*u_filter ;
   deriv_u  =  get_first_derivative(Un,ind(:,4:6),h);
   deriv_u_filter = apply_filter(deriv_u,ind,filter,alpha,mat_alpha);
   M = kappa*kappa* deriv_u_filter.*abs(deriv_u_filter) - apply_filter( deriv_u .*abs(deriv_u) ,ind,filter,alpha,mat_alpha) ;
   csdsq = 0.5 * sum(L.*M) / sum(M.*M); % (Cs * Delta)^2
   dynamic_sub = sqrt(abs(csdsq)) / h ;
endfunction
