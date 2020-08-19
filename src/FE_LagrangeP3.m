function FE_LagrangeP3(N,nu,constant_sub,filter,L,time,nbrpointtemp,Ninterpolation,name,file_spectrum,submethod)
% Solve the 1D forced Burgers equation with cubic Lagrange elements 
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("************************************************************")  
  disp("Finite element cubic Lagrange P3")
  
  length_vec = 3*N;
  h = L/N;% Length of the elements
  DG=4;% Number of Gauss points for numerical integration of the energy
  X = linspace(0,L,3*N)';

  nbrpointtime = nbrpointtemp;
  deltat       = time / nbrpointtime;

  switch submethod
        case 1
             disp("   Full matrix formulation") ;
             Mij = h * [8/105    33/560  -3/140    19/1680  ;
                       33/560   27/70   -27/560  -3/140    ;
                      -3/140   -27/560   27/70    33/560   ;
                       19/1680 -3/140    33/560   8/105     ];
        case 2
             disp("   Lumped matrix formulation") ;
             Mij = h * [210/1680     0         0          0     ;
                         0      630/1680     0          0     ;
                         0         0      630/1680      0     ;
                         0         0         0       210/1680 ];
        otherwise
             disp('Unknown kind of formulation, exiting the code...')
             return
             
  endswitch
  
  MijF = h * [8/105    33/560  -3/140    19/1680  ;
              33/560   27/70   -27/560  -3/140    ;
             -3/140   -27/560   27/70    33/560   ;
              19/1680 -3/140    33/560   8/105     ];
              
  disp("************************************************************")
              
% ******** Rigidity Kij matrix *************  
  Kij = [ 37/10   -189/40   27/20   -13/40  ;
        -189/40   54/5    -297/40   27/20  ;
         27/20   -297/40   54/5    -189/40 ;
        -13/40    27/20   -189/40   37/10        ] * nu / h;
       
  ind = zeros(N,4);
  ind(:,1) = 1:3:3*N-2; 
  ind(:,2) = 2:3:3*N-1; 
  ind(:,3) = 3:3:3*N;
  ind(:,4) = 4:3:3*N+1;
  ind(end,4) = 1;
  
  % Store the indices of the neighbour nodes i-4 i-3 i-2 i-1 i i+1 i+2 i+3 i+4 used in the derivatives
  indfilter = zeros(3*N,9);
  indfilter(:,5) = 1:3*N; % i
  indfilter(:,4) = circshift(indfilter(:,5),1,1) ; % i-1
  indfilter(:,3) = circshift(indfilter(:,5),2,1) ; % i-2
  indfilter(:,2) = circshift(indfilter(:,5),3,1) ; % i-3
  indfilter(:,1) = circshift(indfilter(:,5),4,1) ; % i-4
  indfilter(:,6) = circshift(indfilter(:,5),-1,1); % i+1
  indfilter(:,7) = circshift(indfilter(:,5),-2,1); % i+2
  indfilter(:,8) = circshift(indfilter(:,5),-3,1); % i+3
  indfilter(:,9) = circshift(indfilter(:,5),-4,1); % i+4
  
  dynamic_smag_constant = zeros(nbrpointtemp,1);
  
% ************* Assemble the matrices (with periodic condition) ***********
  M  = zeros(length_vec,length_vec);
  MF = zeros(length_vec,length_vec);
  K  = zeros(length_vec,length_vec);
  for i=1:N
    indmin = 3*i-2;
    indmax = 3*i+1;
    indloc    = 4;
    if (i==(N))
      indmax = 3*i; 
      indloc = 3;
    end
    M(indmin:indmax,indmin:indmax)  += Mij(1:indloc,1:indloc);
    MF(indmin:indmax,indmin:indmax) += MijF(1:indloc,1:indloc);
    K(indmin:indmax,indmin:indmax)  += Kij(1:indloc,1:indloc);
  end
  % last element for periodicity
  M(1,end-2)+= Mij(4,1); M(1,end-1) += Mij(4,2); M(1,end)  += Mij(4,3); M(1,1) += Mij(4,4);
  M(end,1)  += Mij(3,4); M(end-1,1) += Mij(2,4); M(end-2,1)+= Mij(1,4);

  MF(1,end-2)+= MijF(4,1); MF(1,end-1) += MijF(4,2); MF(1,end)  += MijF(4,3); MF(1,1) += MijF(4,4);
  MF(end,1)  += MijF(3,4); MF(end-1,1) += MijF(2,4); MF(end-2,1)+= MijF(1,4);

  K(1,end-2)+= Kij(4,1); K(1,end-1) += Kij(4,2); K(1,end)  += Kij(4,3); K(1,1) += Kij(4,4);
  K(end,1)  += Kij(3,4); K(end-1,1) += Kij(2,4); K(end-2,1)+= Kij(1,4);
  
  M  = sparse(M); % Remove the unnecessary zeros
  MF = sparse(MF);
  K  = sparse(K);
  
% ************* Initialization for numerical integration *****************
% weights for numerical integration (4 points -> 6th order for the subgrid term)
  weight_gauss   = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ];
% coordinates of point for numerical integration (4 points -> 6th order for the subgrid term)
  eta_gauss = [ -0.861136311594953 -0.339981043584856  0.339981043584856  0.861136311594953 ];
  x_gauss = (eta_gauss+1)*h*0.5;  ;
  
% Analytical evaluation of the derivatives of the shape functions at the given integration points
  for i=1:length(weight_gauss)
    d_shape_fct_vector(:,i) = get_deriv_shape_fct(eta_gauss(i)) ;
  end
  
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

  diverged = false;
  for i=2:nbrpointtime+1
%***************** Forcing term with with-noise random phase ******************
    phi2   = 2*pi*rand();    phi3   = 2*pi*rand();
    KForce = 2*pi/L;
    F      = sqrt(deltat) * ( cos(2*KForce*X+phi2) + cos(3*KForce*X+phi3) );
% Uncomment the following line in case of sinus wave with non-linear convection
%      F = 0;
        
%******** Call Runge-Kutta and compute kinematic energy ********
    [u(:,z),dynamic_smag_constant(i-1)] = RK4_FE_Lagrangep3 (u(:,z-1),deltat,h,N,M,MF,K,F,constant_sub,ind,...
                                                             d_shape_fct_vector,weight_gauss,indfilter,filter);
    
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
          xlim([0 1])
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
         
          subplot(2,2,2);
          loglog(0:((length_vec+1)/2-1),spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics,'r','Linewidth',3,...
                 reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3);
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
%  relative_error

  spectralEnergyOut = spectralEnergy(1:((length_vec+1)/2))/nbrPointsStatistics;
  filename2=['Spectral_energy_',name,'.mat'];
  if (~diverged)
    save(filename2,'-ascii','spectralEnergyOut');
  end

%  Uinterpolated = Interpolation(u(:,end),h,N,2,Ninterpolation);%Raffinage de l'affichage
%  tosave = [(0:1/(length(Uinterpolated)):1)' [Uinterpolated;Uinterpolated(1)]];
%  save(strcat(name,'.mat'),'tosave');
        
%  filename=strcat('Energie_',name,'.mat');
%  save(filename,'kinEnergy');
  
end

function [y,smag_sub] = RK4_FE_Lagrangep3 (u,deltat,h,N,M,MF,K,F,constant_sub,ind,...
                                              d_shape_fct_vector,weight,indfilter,filter)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with cubic Lagrange elements
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
  u1 = Un(ind(:,1));  u2 = Un(ind(:,2));    u3 = Un(ind(:,3));    u4 = Un(ind(:,4));
  
%%%%%% Get the Smagorinsky constant in case of dynamic model
  if (filter>0)
     kappa = 2 ; % filter ratio
     smag_sub = get_dynamic_smagorinsky(Un,u1,u2,u3,u4,ind,N,h,kappa,indfilter,filter);
  elseif (filter==0)
     smag_sub = constant_sub ;
  else 
     smag_sub = 0. ;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
% convective term
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid term
  if (smag_sub>0)
    Cj += get_subgrid_terms(smag_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector,weight) ;
  end
  k1 = - M \ (K*Un + Cj);
  Un2 = Un + deltat*0.5*k1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  u1  = Un2(ind(:,1));  u2 = Un2(ind(:,2));    u3 = Un2(ind(:,3));    u4 = Un2(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (smag_sub>0)
    Cj += get_subgrid_terms(smag_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector,weight) ;
  end
  k2 = - M \ (K*Un2 + Cj);
  Un3 = Un + deltat*0.5*k2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  u1  = Un3(ind(:,1));  u2 = Un3(ind(:,2));    u3 = Un3(ind(:,3));    u4 = Un3(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (smag_sub>0)
    Cj += get_subgrid_terms(smag_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector,weight) ;
  end
  k3 = - M \ (K*Un3 + Cj);
  Un4 = Un + deltat*k3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  u1  = Un4(ind(:,1));  u2 = Un4(ind(:,2));    u3 = Un4(ind(:,3));    u4 = Un4(ind(:,4));
% nonlinear terms
  Cj = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N);
% subgrid terms
  if (smag_sub>0)
    Cj += get_subgrid_terms(smag_sub,h,u1,u2,u3,u4,ind,N,d_shape_fct_vector,weight) ;
  end
  k4 = - M \ (K*Un4 + Cj);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 + M\MF*F;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cnonlinear = get_non_linear_lagrange_P3(u1,u2,u3,u4,ind,N)
% Compute the discretization of the nonlinear convection term u * (du/dx)
  Cnonlinear = zeros(3*N,1);
  
  Cnonlinear(ind(:,1)) += (-2240*u1.*u1 + 1611*u1.*u2 + 2025*u2.*u2 - 630*u1.*u3 - 1053*u2.*u3 - 162*u3.*u3 + ...
                           139*u1.*u4 + 180*u2.*u4 - 9*u3.*u4 + 139*u4.*u4)/6720;
                           
  Cnonlinear(ind(:,2)) += (-179*u1.*u1 - 225*u1.*u2 + 72*u1.*u3 + 243*u2.*u3 + 243*u3.*u3 - ...
                           21*u1.*u4 - 18*u2.*u4 - 45*u3.*u4 - 70*u4.*u4)*3/2240 ;
                           
  Cnonlinear(ind(:,3)) += (70*u1.*u1 + 45*u1.*u2 + 18*u1.*u3 - 243*u2.*u2 - 243*u2.*u3 + ...
                           21*u1.*u4 - 72*u2.*u4 + 225*u3.*u4 + 179*u4.*u4)*3/2240;
                           
  Cnonlinear(ind(:,4)) += (2240*u4.*u4 - 1611*u3.*u4 - 2025*u3.*u3 + 630*u2.*u4 + 1053*u2.*u3 + 162*u2.*u2 - ...
                           139*u1.*u4 - 180*u1.*u3 + 9*u1.*u2 - 139*u1.*u1)/6720;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_shape_fct = get_deriv_shape_fct(eta)
% Analytical expression of the derivative of the shape functions in the range -1 < eta < 1
  d_shape_fct = [  1 + 18*eta - 27*eta*eta   ;
                   81*eta*eta - 18*eta - 27  ;
                  -81*eta*eta - 18*eta + 27  ;
                   27*eta*eta + 18*eta - 1   ;
                ] / 16 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_shape_fct = get_deriv_shape_fct_physical(h)
% Analytical expression of the derivative of the shape functions in the range 0 < x < h
  x = linspace(0., h, 4) ;
  eta = 2.*x/h - 1. ;
  d_shape_fct = [  1 + 18*eta - 27*eta.*eta   ;
                   81*eta.*eta - 18*eta - 27  ;
                  -81*eta.*eta - 18*eta + 27  ;
                   27*eta.*eta + 18*eta - 1   ;
                ] / 16 ;
  d_shape_fct = d_shape_fct * 2/h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Csubgrid = get_subgrid_terms(constant_sub,h,u1,u2,u3,u4,ind,N,deriv_shape,weight)
% Compute the subgrid terms needed to model turbulence
  Csubgrid = zeros(3*N,1);
%  factor = constant_sub * constant_sub * 0.125 ;
%  factor = ( constant_sub * h )^2 ;
  factor = constant_sub^2 * 0.5 * h ;
  for i = 1:length(weight) % Loop over integration points
    deriv_u = deriv_shape(1,i) * u1 + deriv_shape(2,i) * u2 + ...
              deriv_shape(3,i) * u3 + deriv_shape(4,i) * u4 ;
              
    factor2 = weight(i) * factor * deriv_u .* abs(deriv_u);
    
    Csubgrid(ind(:,1)) += factor2 * deriv_shape(1,i) ;
    Csubgrid(ind(:,2)) += factor2 * deriv_shape(2,i) ;
    Csubgrid(ind(:,3)) += factor2 * deriv_shape(3,i) ;
    Csubgrid(ind(:,4)) += factor2 * deriv_shape(4,i) ;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smooth = apply_filter(Un,ind,type)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
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
      otherwise
          disp("Unknown type of filter");
          smooth = Un ;
   end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dynamic_sub = get_dynamic_smagorinsky(Un,u1,u2,u3,u4,ind,N,h,kappa,indfilter,filter)
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
   u_filter = apply_filter(Un ,indfilter,filter) ;
   L        = apply_filter(Un.*Un ,indfilter,filter) - u_filter.*u_filter ;
   
  deriv_shape = get_deriv_shape_fct_physical(h) ;
  for i = 1:4
     deriv_u(:,i) = deriv_shape(1,i) * u1 + deriv_shape(2,i) * u2 + ...
                    deriv_shape(3,i) * u3 + deriv_shape(4,i) * u4 ;
  end
  % Average first point of element with last point of previous element
  deriv_u(1:N,1) = 0.5 * ( deriv_u(1:N,1) + deriv_u([N,1:N-1],4) ); 
  
  % Assemble the full array at each node
  deriv_u_full(1:3:3*N-2) = deriv_u(:,1);
  deriv_u_full(2:3:3*N-1) = deriv_u(:,2);
  deriv_u_full(3:3:3*N  ) = deriv_u(:,3);
  
  deriv_u_filter = apply_filter(deriv_u_full,indfilter,filter);
  M = kappa*kappa* deriv_u_filter.*abs(deriv_u_filter) - apply_filter( deriv_u_full .*abs(deriv_u_full) ,indfilter,filter) ;
  
  csdsq = 0.5 *M*L / sum(M.*M) ; % (Cs * Delta)^2
  dynamic_sub = sqrt(abs(csdsq)) / h ;
endfunction