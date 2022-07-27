%--------------------------------------------------------------------------------------------------------
function DGFE_LagrangeP1(N,nu,constant_sub,filter,Alpha_Pade,M_minmod2,L,time,nbrpointtemp,name,file_spectrum)
%--------------------------------------------------------------------------------------------------------
% Solve the 1D forced Burgers equation with cubic Hermite elements
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")
  disp("Discontinuous Galerkin finite element linear Lagrange P1")
%  switch filter
%     case 0
%       disp("   Constant value Smagorinsky model")
%     case 1
%       disp("   Dynamic Smagorinsky model - 3 points stencil for the low-pass filter")
%     case 2
%       disp("   Dynamic Smagorinsky model - 5 points stencil for the low-pass filter")
%     case 3
%       disp("   Dynamic Smagorinsky model - 7 points stencil for the low-pass filter")
%     case 4
%       disp("   Dynamic Smagorinsky model - 9 points stencil for the low-pass filter")
%     case 5
%       disp("   Dynamic Smagorinsky model - Pade low-pass filter")
%     otherwise
       disp("   Direct numerical simulation")
%  end
  disp("************************************************************")

  h=L/N;% Length of the elements
  DG=2;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,N)';
  Xplot(2:2:2*N-2,1) = X(2:end) ; Xplot(3:2:2*N-1,1) = X(2:end) ; Xplot(1,1) = X(1) ;  Xplot(2*N,1) = X(end)+h ;
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;

% ************* Assemble the matrix (with periodic condition) ***********
  Mij = [2  1; 1 2]*h/6;
  M   = zeros(2*N,2*N);
  for i=2:2:2*N-2
    M(i:i+1,i:i+1) += Mij;
  end
  M(1,1) = Mij(2,2); M(1,2*N)= Mij(2,1); M(2*N,2*N) = Mij(1,1); M(2*N,1) = Mij(1,2);
  M  = sparse(M); % Remove the unnecessary zeros

%%%%%%%%%%%%%%%%%%   WARNING %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   LUMP MASS MATRIX USED TO REDUCE THE TIME INCREMENT
%%%%%%%%%%%%%%%%%%   THE DISCRETIZATION IS THEN IDENTICAL TO THE CONSERVATIVE FINITE DIFFERENCE HC2
##  M = eye(2*N) * h * 0.5 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
%  u(1:2:2*N,1) = 2*rand(N,1)-1;
%  u(2:2:2*N,1) = u(1:2:2*N,1) ;
% Sinus solution for non-forced Burgers equation
  u(1:2:2*N,1)=sin(X * 2*pi/L);
  u(2:2:2*N,1) = u(1:2:2*N,1) ;

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
    F(1:2:2*N,1) = sqrt(deltat) * (   cos(2*KForce*X+phi2) +   cos(3*KForce*X+phi3) );
    F(2:2:2*N,1) = F(1:2:2*N,1);
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;

%******** Call Runge-Kutta and compute kinematic energy ********
%    [u(:,z),dynamic_smag_constant(i-1)] = ...
%              RK4_DGFE_Lagrangep1(u(:,z-1),deltat,N,M,nu,h,F,constant_sub,ind,filter,Alpha_Pade ,mat_alpha,M_minmod2);
%    [u(:,z),dynamic_smag_constant(i-1)] = ...
%              RK3SSP_DGFE_Lagrangep1(u(:,z-1),deltat,N,M,nu,h,F,constant_sub,ind,filter,Alpha_Pade ,mat_alpha,M_minmod2);
    [u(:,z),dynamic_smag_constant(i-1)] = ...
              RK5SSP_DGFE_Lagrangep1(u(:,z-1),deltat,N,M,nu,h,F,constant_sub,ind,filter,Alpha_Pade ,mat_alpha,M_minmod2);

    u_avg(1:N,1) = 0.5 * ( u(1:2:2*N-1,z) + u(2:2:2*N,z) );

    kinEnergy(i) = get_kinematic_energy(h,DG,u_avg,N,1);

    if (i*deltat>=timeBeforeStatistics)
      fft_u = fft(u_avg);
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
           filename2=['Spectral_energy_',name,'.mat']; spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end

%        disp(strcat("Max at x=",num2str(X(maxInd)), " and min at x=",num2str(X(minInd))," and distance=",num2str(distance_sinus(i,4))));

        subplot(2,2,1)
        plot(Xplot/L,[uu(2:end) ; uu(1)],'b','Linewidth',3);
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)');
        xlim([0 1]);
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(u_avg)*L/nu)))

%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-u_avg(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
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

%        subplot(2,2,4)
%        mean_Smagorinsky = mean(dynamic_smag_constant(1:i-1));
%        standard_deviation = std(dynamic_smag_constant(1:i-1));
%        standard_deviationp = mean_Smagorinsky + standard_deviation;
%        standard_deviationm = mean_Smagorinsky - standard_deviation;
%        plot((1:(i-1))*deltat,dynamic_smag_constant(1:i-1),'b','Linewidth',3) ; hold on;
%        plot([1 (i-1)]*deltat,[mean_Smagorinsky mean_Smagorinsky],      'r-', 'Linewidth',3);
%        plot([1 (i-1)]*deltat,[standard_deviationp standard_deviationp],'r--','Linewidth',3);
%        plot([1 (i-1)]*deltat,[standard_deviationm standard_deviationm],'r--','Linewidth',3);
%        hold off;
%        grid on; xlabel('Time'); ylabel('Smagorinsky C_s(t)') ;
%        xlim([0 time])

        drawnow;
    end
    z=z+1;

%    CFL=u(:,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(abs(u(:,end)))>1000.)
%    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        diverged = true;
        break;
    end
  end

%  mean_Smagorinsky = mean(dynamic_smag_constant)
%  standard_deviation = std(dynamic_smag_constant)
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

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK3SSP_DGFE_Lagrangep1 (u,deltat,N,M,nu,h,F,constant_sub,ind,filter,alpha,mat_alpha,M_minmod2)
%--------------------------------------------------------------------------------------------------------
% Temporal integration of the 1D Burgers equation with an explicit 3 steps Strong-Stability-Preserving Runge-Kutta scheme
% Spatial discretization with discontinuous linear Lagrange elements
%
% The equation to be solved is
%                  du
%                  -- = f(u,t)
%                  dt
% The explicit 3 steps Strong-Stability-Preserving Runge-Kutta scheme is
% v1     =       U(n) +          deltat * f(U(n))
% v2     = 0.75 *U(n) + 0.25 * ( deltat * f(v1) + v1 )
% U(n+1) = 0.333*U(n) + 0.666* ( deltat * f(v2) + v2 )
%
  nu_over_h = nu / h ;
  smag_sub = 0. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
  Kj = get_viscous_term  (u,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(u,N,ind,smag_sub) ;
  fluxes = get_flux_adv(u,N) + get_flux_viscous(u,N,nu,h);
  k1 = M \ ( Kj - Cj - fluxes);

  Un1 = u + deltat*k1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Kj = get_viscous_term  (Un1,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un1,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un1,N) + get_flux_viscous(Un1,N,nu,h);
  k2 = M \ ( Kj - Cj - fluxes);

  Un2 = 0.75*u + 0.25* ( Un1 + deltat*k2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Kj = get_viscous_term  (Un2,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un2,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un2,N) + get_flux_viscous(Un2,N,nu,h);
  k3 = M \ ( Kj - Cj - fluxes);

  Un3 = u/3.0 + 2.0/3.0 * ( Un2 + deltat*k3 ) ;

  % Apply slope limiter
%  y = Un3 + F ;
  y = slope_lim(Un3,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK5SSP_DGFE_Lagrangep1 (u,deltat,N,M,nu,h,F,constant_sub,ind,filter,alpha,mat_alpha,M_minmod2)
%--------------------------------------------------------------------------------------------------------
% Temporal integration of the 1D Burgers equation with an explicit 5 steps Strong-Stability-Preserving Runge-Kutta scheme
% Spatial discretization with discontinuous linear Lagrange elements
%
  nu_over_h = nu / h ;
  smag_sub = 0. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
  Kj = get_viscous_term  (u,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(u,N,ind,smag_sub) ;
  fluxes = get_flux_adv(u,N) + get_flux_viscous(u,N,nu,h);
  k1 = M \ ( Kj - Cj - fluxes);

  Un1 = u + deltat*0.39175222700392*k1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Kj = get_viscous_term  (Un1,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un1,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un1,N) + get_flux_viscous(Un1,N,nu,h);
  k2 = M \ ( Kj - Cj - fluxes);

  Un2 = 0.44437049406734*u + 0.55562950593266*Un1 + deltat*0.36841059262959*k2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Kj = get_viscous_term  (Un2,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un2,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un2,N) + get_flux_viscous(Un2,N,nu,h);
  k3 = M \ ( Kj - Cj - fluxes);

  Un3 = 0.62010185138540*u + 0.37989814861460*Un2 + deltat*0.25189177424738*k3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  Kj = get_viscous_term  (Un3,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un3,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un3,N) + get_flux_viscous(Un3,N,nu,h);
  k4 = M \ ( Kj - Cj - fluxes);

  Un4 = 0.17807995410773*u + 0.82192004589227*Un3 + deltat*0.54497475021237*k4 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  Kj = get_viscous_term  (Un4,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un4,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un4,N) + get_flux_viscous(Un4,N,nu,h);
  k5 = M \ ( Kj - Cj - fluxes);

  Un5 = 0.00683325884039*u + 0.51723167208978*Un2 + 0.12759831133288*Un3 + 0.34833675773694*Un4 + deltat*0.08460416338212*k4 + deltat*0.22600748319395*k5 ;

  % Apply slope limiter
%  y = Un5 + F ;
  y = slope_lim(Un5,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK4_DGFE_Lagrangep1 (u,deltat,N,M,nu,h,F,constant_sub,ind,filter,alpha,mat_alpha,M_minmod2)
%--------------------------------------------------------------------------------------------------------
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with linear Lagrange elements
%
% The equation to be solved is
%                  du
%                  -- = f(u,t)
%                  dt
% The explicit 4 steps Runge-Kutta scheme is
% U(n+1) = U(n) + 1/6(k1 + 2k2 + 2k3 +k4)
% where k1 = f(U(n)               , t)
%       k2 = f(U(n) + (deltat/2)k1, t + deltat/2)
%       k3 = f(U(n) + (deltat/2)k2, t + deltat/2)
%       k4 = f(U(n) + deltat k3   , t + deltat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
  Un = u;
  nu_over_h = nu / h ;

%%%%%% Get the Smagorinsky constant in case of dynamic model
%  if (filter>0)
%     kappa = 2; % filter ratio
%     smag_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,alpha,mat_alpha);
%  elseif (filter==0)
%     smag_sub = constant_sub ;
%  else
     smag_sub = 0. ;
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  Kj = get_viscous_term  (Un,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un,N) + get_flux_viscous(Un,N,nu,h);
  k1 = M \ ( Kj - Cj - fluxes);

  Un2 = Un + deltat*0.5*k1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  Kj = get_viscous_term  (Un2,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un2,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un2,N) + get_flux_viscous(Un2,N,nu,h);
  k2 = M \ ( Kj - Cj - fluxes);

  Un3 = Un + deltat*0.5*k2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  Kj = get_viscous_term  (Un3,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un3,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un3,N) + get_flux_viscous(Un3,N,nu,h);
  k3 = M \ ( Kj - Cj - fluxes);

  Un4 = Un + deltat*k3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  Kj = get_viscous_term  (Un4,N,ind(:,4:6),nu_over_h) ;
  Cj = get_nonlinear_term(Un4,N,ind,smag_sub) ;
  fluxes = get_flux_adv(Un4,N) + get_flux_viscous(Un4,N,nu,h);
  k4 = M \ ( Kj - Cj - fluxes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 + F;

  y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 ;

  % Apply slope limiter
%  y = y + F ;
  y = slope_lim(y,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function vecC = get_viscous_term(Un,N,ind,nu_over_h)
%--------------------------------------------------------------------------------------------------------
% nu d2(u)/d(x)2
   vec1 = 2:2:2*N-2 ; % indices of first nodes of elements
   vec2 = 3:2:2*N-1 ; % indices of second nodes of elements

   vecC       = zeros(2*N,1) ;
   vecC(1)    = nu_over_h * ( Un(2*N) - Un(1) ) ;     % second node of last element
   vecC(vec1) = nu_over_h * ( Un(vec2) - Un(vec1) ) ; % first node of element
   vecC(vec2) = - vecC(vec1) ; % second node of element
   vecC(2*N)  = - vecC(1)    ; % first node of last element

end

%--------------------------------------------------------------------------------------------------------
function vecC = get_nonlinear_term(Un,N,ind,constant_sub)
%--------------------------------------------------------------------------------------------------------
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
% Convective term - u du/dx
%   vecC = ( Un(ind(:,4)) .* ( Un(ind(:,5)) + Un(ind(:,4)) ) - Un(ind(:,6)) .* ( Un(ind(:,5)) + Un(ind(:,6)) ) ) / 6 ;
   vec1 = 2:2:2*N-2 ; % indices of first nodes of elements
   vec2 = 3:2:2*N-1 ; % indices of second nodes of elements

   vecC       = zeros(2*N,1) ;
   vecC(1)    = ( 2*Un(1)^2     -   Un(2*N)^2   - Un(1)     * Un(2*N)  ) / 6 ; % second node of last element
   vecC(vec1) = (   Un(vec2).^2 - 2*Un(vec1).^2 + Un(vec1) .* Un(vec2) ) / 6 ; % first node of element
   vecC(vec2) = ( 2*Un(vec2).^2 -   Un(vec1).^2 - Un(vec1) .* Un(vec2) ) / 6 ; % second node of element
   vecC(2*N)  = (   Un(1).^2    - 2*Un(2*N).^2  + Un(1)     * Un(2*N)  ) / 6; % first node of last element

% Subgrid term
%   if ( constant_sub>0 )
%      vecC += constant_sub^2 * ( abs(Un(ind(:,6))-Un(ind(:,5))).*(Un(ind(:,6))-Un(ind(:,5))) - ...
%                                 abs(Un(ind(:,5))-Un(ind(:,4))).*(Un(ind(:,5))-Un(ind(:,4))) ) ;
%   endif
end

%--------------------------------------------------------------------------------------------------------
function flux_adv = get_flux_adv(Un,N)
%--------------------------------------------------------------------------------------------------------
   vec1 = 1:2:2*N-1 ; % indices of left values at nodes (minus)
   vec2 = 2:2:2*N   ; % indices of right values at nodes (plus)

   val1 = 0.25 * ( Un(vec2).^2 - Un(vec1).^2 ) ;
   val2 = 0.5 * max( abs([Un(vec1) Un(vec2)]),[],2 ) .* ( Un(vec2) - Un(vec1) ) ;

   % Flux for left values at nodes (minus)
   flux_adv(vec1,1) = val1 - val2 ;

   % Flux for right values at nodes (plus)
   flux_adv(vec2,1) = val1 + val2 ;
end

%--------------------------------------------------------------------------------------------------------
function flux_viscous = get_flux_viscous(Un,N,nu,h)
%--------------------------------------------------------------------------------------------------------
   flux_viscous = zeros(2*N,1) ;

   flux_viscous(3:2:2*N-3) = 0.5 * nu / h * ( Un(1:2:2*N-5) - Un(5:2:2*N-1) ); % Left flux
   flux_viscous(4:2:2*N-2) = 0.5 * nu / h * ( Un(6:2:2*N)   - Un(2:2:2*N-4) ); % Right flux

   % First node
   flux_viscous(1) = 0.5 * nu / h * ( Un(2*N-1) - Un(3)   ); % Left flux
   flux_viscous(2) = 0.5 * nu / h * ( Un(4)     - Un(2*N) ); % Right flux

   % Last node
   flux_viscous(2*N-1) = 0.5 * nu / h * ( Un(2*N-3) - Un(1)     ); % Left flux
   flux_viscous(2*N  ) = 0.5 * nu / h * ( Un(2)     - Un(2*N-2) ); % Right flux
end

%--------------------------------------------------------------------------------------------------------
function u_lim = slope_lim(u,N,h,M_minmod2)
%--------------------------------------------------------------------------------------------------------
  % Precision threshold to apply or not the slope limiter
  % in order to improve accuracy in smooth regions of the solution
  eps0 = 1.0e-6;

  % Average of solution on the elements
  u_mean        = zeros(N,1);
  u_mean(1:N-1) = ( u(2:2:2*N-2) + u(3:2:2*N-1) ) * 0.5;
  u_mean(N)     = ( u(1)         + u(2*N)       ) * 0.5;

  u_tilt        = zeros(N,1);
  u_tilt(1:N-1) = u(3:2:2*N-1) - u_mean(1:N-1);
  u_tilt(N)     = u(1)         - u_mean(N);

  u_tilt2       = u_mean - u(2:2:2*N);

  delta_uplus        = zeros(N,1);
  delta_uplus(1:N-1) = u_mean(2:N) - u_mean(1:N-1);
  delta_uplus(N)     = u_mean(1)   - u_mean(N);

  delta_uminus      = zeros(N,1);
  delta_uminus(2:N) = u_mean(2:N) - u_mean(1:N-1);
  delta_uminus(1)   = u_mean(1)   - u_mean(N);

%  u_tilt_mod  = minmod_vec(u_tilt,  delta_uplus, delta_uminus);
%  u_tilt2_mod = minmod_vec(u_tilt2, delta_uplus, delta_uminus);

  % M is a constant that should be an upper bound on the second derivative at the local extrema
  M = M_minmod2 ;
  u_tilt_mod  = minmod_vec(u_tilt,  delta_uplus + M*h*h*sign(delta_uplus), delta_uminus + M*h*h*sign(delta_uminus));
  u_tilt2_mod = minmod_vec(u_tilt2, delta_uplus + M*h*h*sign(delta_uplus), delta_uminus + M*h*h*sign(delta_uminus));

  u_mod            = zeros(2*N,1);
  u_mod(1)         = u_mean(N)     + u_tilt_mod(1);
  u_mod(3:2:2*N-1) = u_mean(1:N-1) + u_tilt_mod(1:N-1);
  u_mod(2:2:2*N)   = u_mean        - u_tilt2_mod;

  ids = find( abs( u - u_mod ) > eps0 );

  u_lim = u;
  u_lim(ids) = u_mod(ids);
end

%--------------------------------------------------------------------------------------------------------
function r=minmod_vec(a,b,c)
%--------------------------------------------------------------------------------------------------------
   test=[sign(a)==sign(b) sign(b)==sign(c)];
   test2=test(:,1).*test(:,2);
   r=test2.*sign(a).*min([abs(a) abs(b) abs(c)],[],2);
end

%--------------------------------------------------------------------------------------------------------
function dudx = get_first_derivative(Un,ind,h)
%--------------------------------------------------------------------------------------------------------
% ind has size (N,3)
   dudx = ( Un(ind(:,3)) - Un(ind(:,1)) ) * 0.5 / h ;
end

%--------------------------------------------------------------------------------------------------------
function smooth = apply_filter(Un,ind,type,alpha,mat_alpha)
%--------------------------------------------------------------------------------------------------------
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
end

%--------------------------------------------------------------------------------------------------------
function dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,alpha,mat_alpha)
%--------------------------------------------------------------------------------------------------------
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
   u_filter = apply_filter(Un ,ind,filter,alpha,mat_alpha) ;
   L        = apply_filter(Un.*Un ,ind,filter,alpha,mat_alpha) - u_filter.*u_filter ;

   deriv_u  =  get_first_derivative(Un,ind(:,4:6),h);
%   deriv_u_filter = apply_filter(deriv_u,ind,filter,alpha,mat_alpha);
   deriv_u_filter = get_first_derivative(u_filter,ind(:,4:6),h);

   M = kappa*kappa* deriv_u_filter.*abs(deriv_u_filter) - apply_filter( deriv_u .*abs(deriv_u) ,ind,filter,alpha,mat_alpha) ;

   csdsq = 0.5 * sum(L.*M) / sum(M.*M); % (Cs * Delta)^2
   dynamic_sub = sqrt(abs(csdsq)) / h ;
end
