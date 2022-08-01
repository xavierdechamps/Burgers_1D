%--------------------------------------------------------------------------------------------------------
function DGFE_HermiteH3(N,nu,M_DGFE,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum)
%--------------------------------------------------------------------------------------------------------
% Solve the 1D forced Burgers equation with cubic Hermite elements
% The unknowns are the velocity and the first spatial derivative of the velocity, thus 2 unknowns per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")
  disp("Discontinuous finite element cubic Hermite H3")
  disp("   Direct numerical simulation")

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

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(1:4:4*N,1) = 2*rand(N,1)-1; % solution u  left
  u(2:4:4*N,1) = 2*rand(N,1)-1; % solution u' left
  u(3:4:4*N,1) = u(1:4:4*N,1) ; % solution u  right
  u(4:4:4*N,1) = u(2:4:4*N,1) ;  % solution u' right
% Sinus solution for non-forced Burgers equation
%  u(1:4:4*N,1)=sin(X(:) * 2*pi/L) ;
%  u(2:4:4*N,1)=2*pi/L*cos(X(:) * 2*pi/L) ;
%  u(3:4:4*N,1)=u(1:4:4*N,1);
%  u(4:4:4*N,1)=u(2:4:4*N,1);

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
%    u(:,z) = RK4_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,ind,M_DGFE);
%    u(:,z) = RK3SSP_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,ind,M_DGFE);
    u(:,z) = RK5SSP_DGFE_HermiteH3 (u(:,z-1),deltat,h,N,M,MF,K,F,nu,ind,M_DGFE);

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
function y = RK5SSP_DGFE_HermiteH3 (u,deltat,h,N,M,MF,K,F,nu,ind,M_minmod2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First step
% convective term
  Cj = get_non_linear_hermite_H3(u,ind,N,h);
  fluxes = get_flux_adv(u,N) - get_flux_viscous(u,N,nu);

	k1 = - M \ (K*u + Cj + fluxes);
	Un1 = u + deltat*0.39175222700392*k1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% convective term
  Cj = get_non_linear_hermite_H3(Un1,ind,N,h) ;
  fluxes = get_flux_adv(Un1,N) - get_flux_viscous(Un1,N,nu);

	k2 = - M \ (K*Un1 + Cj + fluxes);
	Un2 = 0.44437049406734*u + 0.55562950593266*Un1 + deltat*0.36841059262959*k2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
% convective term
  Cj = get_non_linear_hermite_H3(Un2,ind,N,h) ;
  fluxes = get_flux_adv(Un2,N) - get_flux_viscous(Un2,N,nu);

	k3 = - M \ (K*Un2 + Cj + fluxes);
	Un3 = 0.62010185138540*u + 0.37989814861460*Un2 + deltat*0.25189177424738*k3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
% convective term
  Cj = get_non_linear_hermite_H3(Un3,ind,N,h) ;
  fluxes = get_flux_adv(Un3,N) - get_flux_viscous(Un3,N,nu);

	k4 = - M \ (K*Un3 + Cj + fluxes);
	Un4 = 0.17807995410773*u + 0.82192004589227*Un3 + deltat*0.54497475021237*k4 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
% convective term
  Cj = get_non_linear_hermite_H3(Un4,ind,N,h) ;
  fluxes = get_flux_adv(Un4,N) - get_flux_viscous(Un4,N,nu);

  k5 = - M \ (K*Un3 + Cj + fluxes);
  Un5 = 0.00683325884039*u + 0.51723167208978*Un2 + 0.12759831133288*Un3 + 0.34833675773694*Un4 + deltat*0.08460416338212*k4 + deltat*0.22600748319395*k5 ;

% Apply slope limiter
  y = Un5 + F ;
%  y = slope_lim(Un3,N,h,M_minmod2) + F ;

end

%--------------------------------------------------------------------------------------------------------
function y = RK3SSP_DGFE_HermiteH3 (u,deltat,h,N,M,MF,K,F,nu,ind,M_minmod2)
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
function y = RK4_DGFE_HermiteH3    (u,deltat,h,N,M,MF,K,F,nu,ind,M_minmod2)
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


