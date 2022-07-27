%--------------------------------------------------------------------------------------------------------
function DGFE_HermiteH3(N,nu,constant_sub,Filter,Alpha_Pade,M_DGFE,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
%--------------------------------------------------------------------------------------------------------
% Solve the 1D forced Burgers equation with cubic Hermite elements
% The unknowns are the velocity and the first spatial derivative of the velocity, thus 2 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")
  disp("Finite element cubic Hermite H3")
  switch Filter
     case 0
       disp("   Constant value Smagorinsky model")
##     case 1
##       disp("   Dynamic Smagorinsky model - 3 points stencil for the low-pass filter")
##     case 2
##       disp("   Dynamic Smagorinsky model - 5 points stencil for the low-pass filter")
##     case 3
##       disp("   Dynamic Smagorinsky model - 7 points stencil for the low-pass filter")
##     case 4
##       disp("   Dynamic Smagorinsky model - 9 points stencil for the low-pass filter")
##     case 5
##       disp("   Dynamic Smagorinsky model - Pade low-pass filter")
     otherwise
       disp("   Direct numerical simulation")
  end

  h  = L/N; % Length of the elements
  DG = 4;   % Number of Gauss points for numerical integration of the energy
  X  = zeros(N,1);
  X = linspace(0,L,N)';
  nbrpointtime=nbrpointtemp;
  deltat=time/nbrpointtime;

  Mij=h/35*[  13     11*h/6    9/2   -13*h/12  ;
             11*h/6   h*h/3   13*h/12  -h*h/4  ;
              9/2    13*h/12    13    -11*h/6  ;
            -13*h/12 -h*h/4   -11*h/6   h*h/3 ];

  MijF = Mij ;

  disp("************************************************************")

% ******** Rigidity Kij matrix *************
  Kij = [ 6/5/h    1/10   -6/5/h   1/10  ;
          1/10   2*h/15   -1/10   -h/30  ;
         -6/5/h   -1/10    6/5/h  -1/10  ;
          1/10    -h/30   -1/10  2*h/15] * nu;

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
  M  = zeros(4*N,4*N);
  MF = zeros(4*N,4*N);
  K  = zeros(4*N,4*N);
  for i=3:4:4*N-2
    M(i:i+3,i:i+3)  += Mij;
    MF(i:i+3,i:i+3) += MijF;
    K(i:i+3,i:i+3)  += Kij;
  end
  M(4*N-1:4*N, 4*N-1:4*N) += Mij(1:2,1:2);
  M(4*N-1:4*N, 1:2)       += Mij(1:2,3:4);
  M(1:2      , 1:2)       += Mij(3:4,3:4);
  M(1:2      , 4*N-1:4*N) += Mij(3:4,1:2);

  MF(4*N-1:4*N, 4*N-1:4*N) += MijF(1:2,1:2);
  MF(4*N-1:4*N, 1:2)       += MijF(1:2,3:4);
  MF(1:2      , 1:2)       += MijF(3:4,3:4);
  MF(1:2      , 4*N-1:4*N) += MijF(3:4,1:2);

  K(4*N-1:4*N, 4*N-1:4*N) += Kij(1:2,1:2);
  K(4*N-1:4*N, 1:2)       += Kij(1:2,3:4);
  K(1:2      , 1:2)       += Kij(3:4,3:4);
  K(1:2      , 4*N-1:4*N) += Kij(3:4,1:2);

  M  = sparse(M); % Remove the unnecessary zeros
  MF = sparse(MF);
  K  = sparse(K);

% ************* Initialization for numerical integration *****************
% weights for numerical integration (5 points -> 8th order for the subgrid term)
%  weight_gauss   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ]; % 4 points
  weight_gauss   = [ 0.2369268850561891  0.4786286704993665  0.5688888888888889  0.4786286704993665  0.2369268850561891  ]; % 5 points
% coordinates of point for numerical integration (5 points -> 8th order for the subgrid term)
%  eta_gauss = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ]; % 4 points
  eta_gauss = [ -0.9061798459386640  -0.5384693101056831  0.0000000000000000  0.5384693101056831  0.9061798459386640 ]; % 5 points
  x_gauss = (eta_gauss+1)*h*0.5;  ;

% Analytical evaluation of the derivatives of the shape functions at the given integration points
  for i=1:length(weight_gauss)
    d_shape_fct_vector(:,i) = get_deriv_shape_fct(eta_gauss(i)) ;
  end

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
%  u(1:4:4*N,1) = 2*rand(N,1)-1; % solution u  left
%  u(2:4:4*N,1) = 2*rand(N,1)-1; % solution u' left
%  u(3:4:4*N,1) = u(1:4:4*N,1) ; % solution u  right
%  u(4:4:4*N,1) = u(2:4:4*N,1) ;  % solution u' right
% Sinus solution for non-forced Burgers equation
  u(1:4:4*N,1)=sin(X(:) * 2*pi/L) ;
  u(2:4:4*N,1)=2*pi/L*cos(X(:) * 2*pi/L) ;
  u(3:4:4*N,1)=u(1:4:4*N,1);
  u(4:4:4*N,1)=u(2:4:4*N,1);

% Store the indices of the neighbour nodes i-4 i-3 i-2 i-1 i i+1 i+2 i+3 i+4 used in the derivatives
  indfilter = zeros(N,9);
  indfilter(:,5) = 1:N; % i
  indfilter(:,4) = circshift(indfilter(:,5),1,1) ; % i-1
  indfilter(:,3) = circshift(indfilter(:,5),2,1) ; % i-2
  indfilter(:,2) = circshift(indfilter(:,5),3,1) ; % i-3
  indfilter(:,1) = circshift(indfilter(:,5),4,1) ; % i-4
  indfilter(:,6) = circshift(indfilter(:,5),-1,1); % i+1
  indfilter(:,7) = circshift(indfilter(:,5),-2,1); % i+2
  indfilter(:,8) = circshift(indfilter(:,5),-3,1); % i+3
  indfilter(:,9) = circshift(indfilter(:,5),-4,1); % i+4
  dynamic_smag_constant = zeros(nbrpointtemp,1);
  mat_alpha = zeros(N,N) ;
  for i=1:N
     mat_alpha(i, indfilter(i,4:6)) = [Alpha_Pade , 1 , Alpha_Pade] ;
  end
  mat_alpha = sparse(mat_alpha);

  timeBeforeStatistics = 10;

  j=2;  z=2;  ind_error=1;
  filename=[name,num2str(1),'.mat'];
  uu=u(:,1);
%save(filename,'uu');

  kinEnergy = zeros(1:nbrpointtime+1,1); kinEnergy(1) = get_kinematic_energy(h,DG,u(:,1),N,2); kinEnergyMean = 0; nbrPointsStatistics = 0;
  spectralEnergy=zeros(Ninterpolation*N,1);
  reference_spectrum=load(file_spectrum);
  diverged = false;
  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2 = 2*pi*rand();    phi3 = 2*pi*rand();

    KForce = 2*pi/L;
    F(1:4:4*N,1)      =   sqrt(deltat) *         (   cos(2*KForce*X+phi2) +   cos(3*KForce*X+phi3) );
    F(2:4:4*N,1)      = - sqrt(deltat) * KForce* ( 2*sin(2*KForce*X+phi2) + 3*sin(3*KForce*X+phi3) );
    F(3:4:4*N,1)      = F(1:4:4*N,1);
    F(4:4:4*N,1)      = F(2:4:4*N,1);
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;

%******** Call Runge-Kutta and compute kinematic energy ********
%    [u(:,z),dynamic_smag_constant(i-1)] = ...
%             RK4_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,constant_sub,Filter,indfilter,...
%                               Alpha_Pade ,mat_alpha,ind,d_shape_fct_vector,weight_gauss,M_DGFE);
%    [u(:,z),dynamic_smag_constant(i-1)] = ...
%             RK3SSP_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,constant_sub,Filter,indfilter,...
%                               Alpha_Pade ,mat_alpha,ind,d_shape_fct_vector,weight_gauss,M_DGFE);
    [u(:,z),dynamic_smag_constant(i-1)] = ...
             RK5SSP_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,constant_sub,Filter,indfilter,...
                               Alpha_Pade ,mat_alpha,ind,d_shape_fct_vector,weight_gauss,M_DGFE);

    u_avg(1:2:2*N,1) = 0.5 * ( u(1:4:4*N,z) + u(3:4:4*N,z) );
    u_avg(2:2:2*N,1) = 0.5 * ( u(2:4:4*N,z) + u(4:4:4*N,z) );

    kinEnergy(i) = get_kinematic_energy(h,DG,u_avg,N,2);

    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1) ;
      Uinterpolated        = Interpolation(u_avg,h,N,1,Ninterpolation); % Interpolate within each element to show internal oscillations
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

        Uinterpolated = Interpolation(u_avg,h,N,1,Ninterpolation); % Interpolate within each element to show internal oscillations

        subplot(2,2,1)
        plot(0:1/(length(Uinterpolated)):1,[Uinterpolated;Uinterpolated(1)],'b','Linewidth',3)
        title(strcat('Time= ',num2str(i*deltat),', Re= ',num2str(mean(Uinterpolated)*L/nu)))
        xlabel('x/(2*\pi)'); ylabel('u(t)'); grid on
        xlim([0 1])

%        hold on
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,[0:(2/(N*Ninterpolation)):1],100)' ;
%        plot([0:(2/(N*Ninterpolation)):1]/2,sol_theory,'r')
%        hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-Uinterpolated(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;

%        hold on
%        plot([X;X(end)+h]/L,[u(1:2:2*N-1,end);u(1,end)],'b*')
%        hold off

        subplot(2,2,3);
        plot((1:i)*deltat,kinEnergy(1:i),'b','Linewidth',3);
        grid on; xlabel('Time'); ylabel('E(t)');
        xlim([0 time])

        subplot(2,2,2);
        loglog(0:N-1,spectralEnergy(1:N)/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
        grid on; xlabel('k'); ylabel('E(k)');
        xlim([1 reference_spectrum(end,1)])

##        subplot(2,2,4)
##        mean_Smagorinsky = mean(dynamic_smag_constant(1:i-1));
##        standard_deviation = std(dynamic_smag_constant(1:i-1));
##        standard_deviationp = mean_Smagorinsky + standard_deviation;
##        standard_deviationm = mean_Smagorinsky - standard_deviation;
##        plot((1:(i-1))*deltat,dynamic_smag_constant(1:i-1),'b','Linewidth',3) ; hold on;
##        plot([1 (i-1)]*deltat,[mean_Smagorinsky mean_Smagorinsky],      'r-', 'Linewidth',3);
##        plot([1 (i-1)]*deltat,[standard_deviationp standard_deviationp],'r--','Linewidth',3);
##        plot([1 (i-1)]*deltat,[standard_deviationm standard_deviationm],'r--','Linewidth',3);
##        hold off;
##        grid on; xlabel('Time'); ylabel('Smagorinsky C_s(t)') ;
##        xlim([0 time])

        drawnow ;
    end

    z=z+1;
%    CFL=u(1:2:end,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(abs(u(1:2:end,end)))>1000.)
%    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        diverged = true;
        break;
    end
  end

%  relative_error

  spectralEnergyOut = spectralEnergy(1:N)/nbrPointsStatistics;
%  filename2=['Spectral_energy_',name,'.mat'];
  filename2=['Spectral_energy_',name,'.mat'];
  if (~diverged)
    save(filename2,'-ascii','spectralEnergyOut');
  end

%  Uinterpolated = Interpolation(u(:,end),h,N,1,Ninterpolation);
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');

%  filename=strcat('Energy_',name,'.mat');
%  save(filename,'kinEnergy');

end

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK5SSP_DGFE_HermiteH3 (u,deltat,h,N,M,MF,K,F,nu,constant_sub,filter,indfilter,...
                                            alpha,mat_alpha,ind,d_shape_fct_vector,weight,M_minmod2)
%--------------------------------------------------------------------------------------------------------
    if (filter==0)
     smag_sub = constant_sub ;
  else
     smag_sub = 0. ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
% convective term
  Cj = get_non_linear_hermite_H3(u,ind,N,h);
  fluxes = get_flux_adv(u,N) - get_flux_viscous(u,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,u,ind,N,d_shape_fct_vector,weight);
%  end

	k1 = - M \ (K*u + Cj + fluxes);
	Un1 = u + deltat*0.39175222700392*k1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% convective term
  Cj = get_non_linear_hermite_H3(Un1,ind,N,h) ;
  fluxes = get_flux_adv(Un1,N) - get_flux_viscous(Un1,N,nu);
% subgrid term
%  if (smag_sub>0)
%    Cj += get_subgrid_terms(smag_sub,h,Un1,ind,N,d_shape_fct_vector,weight);
%  end

	k2 = - M \ (K*Un1 + Cj + fluxes);
	Un2 = 0.44437049406734*u + 0.55562950593266*Un1 + deltat*0.36841059262959*k2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
% convective term
  Cj = get_non_linear_hermite_H3(Un2,ind,N,h) ;
  fluxes = get_flux_adv(Un2,N) - get_flux_viscous(Un2,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un2,ind,N,d_shape_fct_vector,weight);
%	end

	k3 = - M \ (K*Un2 + Cj + fluxes);
	Un3 = 0.62010185138540*u + 0.37989814861460*Un2 + deltat*0.25189177424738*k3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
% convective term
  Cj = get_non_linear_hermite_H3(Un3,ind,N,h) ;
  fluxes = get_flux_adv(Un3,N) - get_flux_viscous(Un3,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un3,ind,N,d_shape_fct_vector,weight);
%	end

	k4 = - M \ (K*Un3 + Cj + fluxes);
	Un4 = 0.17807995410773*u + 0.82192004589227*Un3 + deltat*0.54497475021237*k4 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
% convective term
  Cj = get_non_linear_hermite_H3(Un4,ind,N,h) ;
  fluxes = get_flux_adv(Un4,N) - get_flux_viscous(Un4,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un4,ind,N,d_shape_fct_vector,weight);
%	end

	k5 = - M \ (K*Un3 + Cj + fluxes);
	Un5 = 0.00683325884039*u + 0.51723167208978*Un2 + 0.12759831133288*Un3 + 0.34833675773694*Un4 + deltat*0.08460416338212*k4 + deltat*0.22600748319395*k5 ;

% Apply slope limiter
  y = Un5 + F ;
%  y = slope_lim(Un3,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK3SSP_DGFE_HermiteH3 (u,deltat,h,N,M,MF,K,F,nu,constant_sub,filter,indfilter,...
                                            alpha,mat_alpha,ind,d_shape_fct_vector,weight,M_minmod2)
%--------------------------------------------------------------------------------------------------------
  if (filter==0)
     smag_sub = constant_sub ;
  else
     smag_sub = 0. ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
% convective term
  Cj = get_non_linear_hermite_H3(u,ind,N,h);
  fluxes = get_flux_adv(u,N) - get_flux_viscous(u,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,u,ind,N,d_shape_fct_vector,weight);
%  end

	k1 = - M \ (K*u + Cj + fluxes);
	Un1 = u + deltat*k1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% convective term
  Cj = get_non_linear_hermite_H3(Un1,ind,N,h) ;
  fluxes = get_flux_adv(Un1,N) - get_flux_viscous(Un1,N,nu);
% subgrid term
%  if (smag_sub>0)
%    Cj += get_subgrid_terms(smag_sub,h,Un1,ind,N,d_shape_fct_vector,weight);
%  end

	k2 = - M \ (K*Un1 + Cj + fluxes);
	Un2 = 0.75*u + 0.25* ( Un1 + deltat*k2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
% convective term
  Cj = get_non_linear_hermite_H3(Un2,ind,N,h) ;
  fluxes = get_flux_adv(Un2,N) - get_flux_viscous(Un2,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un2,ind,N,d_shape_fct_vector,weight);
%	end

	k3 = - M \ (K*Un2 + Cj + fluxes);
	Un3 = u/3.0 + 2.0/3.0 * ( Un2 + deltat*k3 ) ;

% Apply slope limiter
  y = Un3 + F ;
%  y = slope_lim(Un3,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function [y,smag_sub] = RK4_DGFE_HermiteH3 (u,deltat,h,N,M,MF,K,F,nu,constant_sub,filter,indfilter,...
                                            alpha,mat_alpha,ind,d_shape_fct_vector,weight,M_minmod2)
%--------------------------------------------------------------------------------------------------------
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

%%%%%% Get the Smagorinsky constant in case of dynamic model
##  if (filter>0)
##     kappa = 2; % filter ratio
##     smag_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,indfilter,alpha,mat_alpha);
##  elseif (filter==0)
  if (filter==0)
     smag_sub = constant_sub ;
  else
     smag_sub = 0. ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% convective term
  Cj = get_non_linear_hermite_H3(Un,ind,N,h);
  fluxes = get_flux_adv(Un,N) - get_flux_viscous(Un,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un,ind,N,d_shape_fct_vector,weight);
%  end

	k1 = - M \ (K*Un + Cj + fluxes);
	Un2 = Un + deltat*0.5*k1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
% convective term
  Cj = get_non_linear_hermite_H3(Un2,ind,N,h) ;
  fluxes = get_flux_adv(Un2,N) - get_flux_viscous(Un2,N,nu);
% subgrid term
%  if (smag_sub>0)
%    Cj += get_subgrid_terms(smag_sub,h,Un2,ind,N,d_shape_fct_vector,weight);
%  end

	k2 = - M \ (K*Un2 + Cj + fluxes);
	Un3 = Un + deltat*0.5*k2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
% convective term
  Cj = get_non_linear_hermite_H3(Un3,ind,N,h) ;
  fluxes = get_flux_adv(Un3,N) - get_flux_viscous(Un3,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un3,ind,N,d_shape_fct_vector,weight);
%	end

	k3 = - M \ (K*Un3 + Cj + fluxes);
	Un4 = Un + deltat*k3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
% convective term
  Cj = get_non_linear_hermite_H3(Un4,ind,N,h) ;
  fluxes = get_flux_adv(Un4,N) - get_flux_viscous(Un4,N,nu);
% subgrid term
%  if (smag_sub>0)
%     Cj += get_subgrid_terms(smag_sub,h,Un4,ind,N,d_shape_fct_vector,weight);
% end

	k4 = - M \ (K*Un4 + Cj + fluxes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 + M\MF*F;
	y = Un + deltat*( k1 + 2*k2 + 2*k3 + k4 )/6 ;

% Apply slope limiter
  y = y + F ;
%  y = slope_lim(y,N,h,M_minmod2) + F ;

end


%--------------------------------------------------------------------------------------------------------
function Cnonlinear = get_non_linear_hermite_H3(u,ind,N,h)
%--------------------------------------------------------------------------------------------------------
% Compute the discretization of the nonlinear convection term u * (du/dx)
	Cnonlinear = zeros(4*N,1);

  vec1  = 3:4:4*N-1;               % Indices of first nodes in elements
  vec1d = 4:4:4*N;                 % First nodes of elements, derivatives
  vec2  = 5:4:4*N; vec2(end+1)=1;  % Indices of second nodes in elements
  vec2d = 6:4:4*N; vec2d(end+1)=2; % Second nodes of elements, derivatives

  u1  = u(vec1);      u2  = u(vec2);
  du1 = u(vec1d)*h;   du2 = u(vec2d)*h;

  Cnonlinear(vec1) = ( u1.*(-280*u1 + 50*du1 + 140*u2 - 34*du2) + du1.*(5*du1 + 34*u2 - 8*du2 ) + u2.*(140*u2 - 50*du2 )  + 5*du2.^2 )/840;
  Cnonlinear(vec1d)= ( u1.*(50*u1 + 5*du1 - 16*u2 + 3*du2)      + du1.*(-5*u2 + du2 )           + u2.*(-34*u2 + 11*du2 )  -   du2.^2 )*(-h/840);
  Cnonlinear(vec2) = ( u1.*(140*u1 + 50*du1 + 140*u2 - 34*du2)  + du1.*(5*du1 + 34*u2 - 8*du2 ) + u2.*(-280*u2 - 50*du2 ) + 5*du2.^2 )*(-1/840);
  Cnonlinear(vec2d)= ( u1.*(34*u1 + 11*du1 + 16*u2 - 5*du2)     + du1.*(du1 + 3*u2 - du2 )      + u2.*(-50*u2 + 5*du2 )              )*(h/840);
end

%--------------------------------------------------------------------------------------------------------
function flux_adv = get_flux_adv(Un,N)
%--------------------------------------------------------------------------------------------------------
   vec1 = 1:4:4*N ; % indices of left values at nodes (minus)
   vec2 = 3:4:4*N ; % indices of right values at nodes (plus)

   val1 = 0.25 * ( Un(vec2).^2 - Un(vec1).^2 ) ;
   val2 = 0.5 * max( abs([Un(vec1) Un(vec2)]),[],2 ) .* ( Un(vec2) - Un(vec1) ) ;

   flux_adv = zeros(4*N,1) ;

   % Flux for left values at nodes (minus)
   flux_adv(vec1) = val1 - val2 ;

   % Flux for right values at nodes (plus)
   flux_adv(vec2) = val1 + val2 ;
end

%--------------------------------------------------------------------------------------------------------
function flux_viscous = get_flux_viscous(Un,N,nu)
%--------------------------------------------------------------------------------------------------------
   flux_viscous = zeros(4*N,1) ;

   flux_viscous(1:4:4*N) = 0.5 * nu * ( Un(2:4:4*N) + Un(4:4:4*N) ); % Left flux for u
   flux_viscous(2:4:4*N) = 0.5 * nu * ( Un(3:4:4*N) - Un(1:4:4*N) ); % Left flux for u'
   flux_viscous(3:4:4*N) = - flux_viscous(1:4:4*N)                 ; % Right flux for u
   flux_viscous(4:4:4*N) =   flux_viscous(2:4:4*N)                 ; % Right flux for u'
end

%--------------------------------------------------------------------------------------------------------
function u_lim = slope_lim(u,N,h,M_minmod2)
%--------------------------------------------------------------------------------------------------------
  % Precision threshold to apply or not the slope limiter
  % in order to improve accuracy in smooth regions of the solution
  eps0 = 1.0e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLOPE LIMITING THE UNKNOWN u
  % Average of solution on the elements
  u_mean        = zeros(N,1);
  u_mean(1:N-1) = ( u(3:4:4*N-5) + u(5:4:4*N-3) ) * 0.5;
  u_mean(N)     = ( u(1)         + u(4*N-1)     ) * 0.5;

  u_tilt        = zeros(N,1);
  u_tilt(1:N-1) = u(5:4:4*N-3) - u_mean(1:N-1);
  u_tilt(N)     = u(1)         - u_mean(N);

  u_tilt2       = u_mean - u(3:4:4*N-1);

  delta_uplus        = zeros(N,1);
  delta_uplus(1:N-1) = u_mean(2:N) - u_mean(1:N-1);
  delta_uplus(N)     = u_mean(1)   - u_mean(N);

  delta_uminus      = zeros(N,1);
  delta_uminus(2:N) = u_mean(2:N) - u_mean(1:N-1);
  delta_uminus(1)   = u_mean(1)   - u_mean(N);

  % M is a constant that should be an upper bound on the second derivative at the local extrema
  M = M_minmod2 ;
  u_tilt_mod  = minmod_vec(u_tilt,  delta_uplus + M*h*h*sign(delta_uplus), delta_uminus + M*h*h*sign(delta_uminus));
  u_tilt2_mod = minmod_vec(u_tilt2, delta_uplus + M*h*h*sign(delta_uplus), delta_uminus + M*h*h*sign(delta_uminus));

  u_mod            = zeros(2*N,1);
  u_mod(1)         = u_mean(N)     + u_tilt_mod(1);
  u_mod(3:2:2*N-1) = u_mean(1:N-1) + u_tilt_mod(1:N-1);
  u_mod(2:2:2*N)   = u_mean        - u_tilt2_mod;

  ids = find( abs( u(1:2:4*N) - u_mod ) > eps0 ) ;

  u_lim = u;
  u_lim(2*ids-1) = u_mod(ids);

  if (~isempty(ids))
    left_nodes = ids ( find( mod(ids,2)==0 & ids~=2*N ) ) ; % exclude last node
    right_nodes = ids ( find( mod(ids,2)==1  & ids~=1 ) ) ; % exclude first node
    test = right_nodes(1) == 1 ; % first node
    test2= left_nodes(end) == 2*N; % last node

    u_lim(2*left_nodes) = ( u_lim(2*left_nodes+1) - u_lim(2*left_nodes-1) )/h ; % limit slope for u+'
    u_lim(2*right_nodes) = ( u_lim(2*right_nodes-1) - u_lim(2*right_nodes-3) )/h ; % limit slope for u-'
    if (test)
      u_lim(2) = ( u_lim(1) - u_lim(4*N-1) )/h;
    endif
    if (test2)
      u_lim(4*N) = ( u_lim(1) - u_lim(4*N-1) )/h;
    endif
  endif

end

%--------------------------------------------------------------------------------------------------------
function r=minmod_vec(a,b,c)
%--------------------------------------------------------------------------------------------------------
   test=[sign(a)==sign(b) sign(b)==sign(c)];
   test2=test(:,1).*test(:,2);
   r=test2.*sign(a).*min([abs(a) abs(b) abs(c)],[],2);
end

%--------------------------------------------------------------------------------------------------------
function d_shape_fct = get_deriv_shape_fct(ksi)
%--------------------------------------------------------------------------------------------------------
% Analytical expression of the derivative of the shape functions
	d_shape_fct = [ 0.5 * (ksi - 2)*(1 + ksi) + 0.25 * (1 + ksi)^2                             ;
	                0.5 * (1 + ksi)*((1 + ksi)*0.5 - 1) + 0.5* (0.25* (1 + ksi)^2 -ksi)        ;
	                0.5 * (2 - ksi)*(1 + ksi) - 0.25* (1 + ksi)^2                              ;
	                0.25* (ksi -1)*(1 + ksi) + 0.125* (1 + ksi)^2
	              ];
end

%--------------------------------------------------------------------------------------------------------
function Csubgrid = get_subgrid_terms(constant_sub,h,Un,ind,N,deriv_shape,weight)
%--------------------------------------------------------------------------------------------------------
% Compute the Smagorinsky subgrid term

	Csubgrid = zeros(2*N,1);
%	u1 = Un(ind(:,1));	du1 = Un(ind(:,2)); % not needed for the numerical integration
  u2 = Un(ind(:,3));  du2 = Un(ind(:,4));
  u3 = Un(ind(:,5));  du3 = Un(ind(:,6));

  factor = constant_sub^2 * 2 * h ;
%  factor = constant_sub^2 * 4 ;
  factor_deriv = 1 ; %  h * 0.5 ;

	for i = 1:length(weight) % Loop over the integration points
		deriv_u = deriv_shape(1,i) * u2 + deriv_shape(2,i) * du2 * factor_deriv  + ...
              deriv_shape(3,i) * u3 + deriv_shape(4,i) * du3 * factor_deriv  ;

	  factor2 = weight(i) * factor * deriv_u .* abs(deriv_u);

		Csubgrid(ind(:,3)) += factor2 * deriv_shape(1,i) ;
%		Csubgrid(ind(:,4)) += factor2 * deriv_shape(2,i) ;
		Csubgrid(ind(:,5)) += factor2 * deriv_shape(3,i) ;
%		Csubgrid(ind(:,6)) += factor2 * deriv_shape(4,i) ;
	end

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
function dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,indfilter,alpha,mat_alpha)
%--------------------------------------------------------------------------------------------------------
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
   u_filter = apply_filter(Un(1:2:end) ,indfilter,filter,alpha,mat_alpha) ;
   L        = apply_filter(Un(1:2:end).*Un(1:2:end) ,indfilter,filter,alpha,mat_alpha) - u_filter.*u_filter ;

   deriv_u_filter = apply_filter(Un(2:2:end),indfilter,filter,alpha,mat_alpha);
   M = kappa*kappa* deriv_u_filter.*abs(deriv_u_filter) - apply_filter( Un(2:2:end) .*abs(Un(2:2:end)) ,indfilter,filter,alpha,mat_alpha) ;

   csdsq = 0.5 * sum(L.*M) / sum(M.*M); % (Cs * Delta)^2
   dynamic_sub = sqrt(abs(csdsq)) / h ;
endfunction
