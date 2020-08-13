function FD_compact_spectral (N,nu,constant_sub,filter,L,time,nbrpointtemp,name,file_spectrum,submethod)
% Solve the 1D forced Burgers equation with compact finite difference schemes
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp('*********************************************************************')  
  disp('Finite difference - Compact schemes with spectral-like resolution')

  h = L/(N);% Length of the elements
  X = linspace(0,L,N)';
  nbrpointtime = nbrpointtemp;
  deltat       = time / nbrpointtime;

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1) = 2*rand(N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(:,1) = sin(X * 2*pi/L);

  timeBeforeStatistics = 10;

% Multiplicative factors for the first and second spatial derivatives
% The first column contains the constants for the first spatial derivative
% The second column contains the constants for the second spatial derivative
% The constants in the rows are respectively [ alpha  beta  a  b  c ]'
  factors_deriv = zeros(5,2) ;
  switch submethod
     case 1
% Optimal - see the doc
        factors_deriv(1:5,1) = [0.5771439  0.0896406  1.3025166  0.99355  0.03750245]';
        factors_deriv(1:5,2) = [0.50209266 0.05569169 0.21564935 1.723322 0.17659730]';
        disp('   Optimal order')
     case 2
% Order 10 - see the doc
        factors_deriv(1:5,1) = [0.5  1/20  17/12  101/150  0.01]';
        factors_deriv(1:5,2) = [334/899 43/1798 1065/1798 1038/899 79/1798]';
        disp('   10th order')
     case 3
% Order 8 - see the doc
        factors_deriv(1:5,1) = [4/9 1/36 40/27 25/54 0]';
        alpha = 344/1179;
        factors_deriv(1:5,2) = [alpha (38*alpha-9)/214 (696-1191*alpha)/428 (2454*alpha-294)/535 0]';
        disp('   8th order')
     case 4
% Order 6 - see the doc
        alpha = 1/3;  beta=0;
        factors_deriv(1:5,1) = [alpha beta (9+alpha-20*beta)/6 (-9+32*alpha+62*beta)/15 (1-3*alpha+12*beta)/10]';
        alpha = 2/11; beta=0;
        factors_deriv(1:5,2) = [alpha beta (6-9*alpha-12*beta)*0.25 (-3+24*alpha-6*beta)*0.2 (2-11*alpha+124*beta)/20]';
        disp('   6th order')
     case 5
% Order 4 central scheme - see the doc
        alpha = 0; beta=0;
        factors_deriv(1:5,1) = [0 0 2*(2+alpha)/3 (4*alpha-1)/3  0]';
        alpha = 0; beta=0;
        factors_deriv(1:5,2) = [0 0 4*(1-alpha)/3 (10*alpha-1)/3  0]';
        disp('   4th order')
     case 6
% Order 2 central scheme - see the doc
        factors_deriv(1:5,1) = [0 0 1 0  0]';
        factors_deriv(1:5,2) = [0 0 1 0  0]';
        disp('   2nd order')
     otherwise
        disp('Unknown kind of compact scheme, exiting the code...')
        return
  end
  disp('*********************************************************************')  
  
% Scale the factors
  factors_deriv(3,1) = factors_deriv(3,1)/(2*h); factors_deriv(4,1) = factors_deriv(4,1)/(4*h);   factors_deriv(5,1) = factors_deriv(5,1)/(6*h); 
  factors_deriv(3,2) = factors_deriv(3,2)/(h*h); factors_deriv(4,2) = factors_deriv(4,2)/(4*h*h); factors_deriv(5,2) = factors_deriv(5,2)/(9*h*h); 
  factors_deriv(3:5,2) = factors_deriv(3:5,2) * nu ;
  
% Store the indices of the neighbour nodes i-3 i-2 i-1 i i+1 i+2 i+3 used in the derivatives
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
  
% Create the matrices to solve the system of equations for the first 
% and second spatial derivatives
   mat_deriv1=zeros(N,N) ;
   mat_deriv2=zeros(N,N) ;
   for i=1:N
      mat_deriv1(i, ind(i,3:7)) = [factors_deriv(2,1) , factors_deriv(1,1) , 1 , factors_deriv(1,1) , factors_deriv(2,1)] ;
      mat_deriv2(i, ind(i,3:7)) = [factors_deriv(2,2) , factors_deriv(1,2) , 1 , factors_deriv(1,2) , factors_deriv(2,2)] ;
   endfor
   mat_deriv1 = sparse(mat_deriv1);
   mat_deriv2 = sparse(mat_deriv2);
      
  kinEnergy = zeros(nbrpointtemp+1,1); kinEnergy(1) = u(:,1)' * u(:,1) * h * 0.5 ;
% energy_conv is the numerical energy produced by the spatial discretization of the convective term,
% it must be equal to zero for an energy conservative scheme
  energy_conv = zeros(nbrpointtemp+1,1);
    
%  filename=strcat(name,num2str(1),'.mat'); uu=u(:,1); save(filename,'uu');

  z=2; j=2; ind_error=1;
  nbrPointsStatistics=0;  kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);

  dynamic_smag_constant = zeros(nbrpointtemp,1);
  diverged = false;
  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand();    phi3=2*pi*rand();
    KForce=2*pi/L;
    F = sqrt(deltat)*(cos(2*KForce*X+phi2)+cos(3*KForce*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;
    
%******** Call Runge-Kutta and compute kinematic energy ********
    [u(:,z), energy_conv(i),dynamic_smag_constant(i-1)] = ...
              RK4_FD_compact_spectral (u(:,z-1),deltat,N,mat_deriv1,mat_deriv2,F,h,constant_sub,factors_deriv,ind,filter);
    kinEnergy(i) = h*0.5*u(:,z)'*u(:,z);
        
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);
      fft_u = fft(u(:,end));
      newSpectralEnergy = fft_u.*conj(fft_u)/N/N;
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
        
        reynolds_number        = mean(uu)*L/nu;
        kolmogorov_length      = L*reynolds_number^(-0.75);
        kolmogorov_wave_length = 2*pi/kolmogorov_length ;
                
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ));
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'-ascii','spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ));
        end
        
        subplot(2,2,1);
        plot(X/L, uu,'b','Linewidth',3)
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        xlim([0 1])
        title(strcat('Time= ',num2str((i-1)*deltat),', Re= ',num2str(mean(uu)*L/nu)))
        
%        sol_theory = get_analytical_solution((i-1)*deltat,nu,X(1:(end/2 + 1)),100) ;
%        hold on; plot(X(1:(end/2 + 1))/L,sol_theory,'r'); hold off
%        relative_error(ind_error,1) = (i-1)*deltat;
%        relative_error(ind_error,2) = sqrt( sum( (sol_theory-uu(1:length(sol_theory))).^2 ) ) / sqrt( sum( sol_theory.^2 ) ) ;
%        disp(strcat('Time : ',num2str(relative_error(ind_error,1)),' and relative error : ', num2str(relative_error(ind_error,2)) ))
%        ind_error = ind_error + 1;
        
        subplot(2,2,3)
        plot((0:(i-1))*deltat,kinEnergy(1:i),'b','Linewidth',3)
        grid on; xlabel('Time'); ylabel('E(t)')
        xlim([0 time])
        
        subplot(2,2,2)
%        subplot(1,2,2)
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
        grid on; xlabel('k'); ylabel('E(k)')
        xlim([1 reference_spectrum(end,1)])
        
%        subplot(2,2,4)
%        plot((0:(i-1))*deltat,energy_conv(1:i),'b','Linewidth',3)
%        grid on; xlabel('Time'); ylabel('E_{prod}(t)')
%        xlim([0 time])
        
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
  
    spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
    filename2=strcat('Spectral_energy_',name,'.mat');
    if (~diverged)
      save(filename2,'-ascii','spectralEnergyOut');
    end
  %filename=strcat('Energy_',name,'.mat');
  %save(filename,'kinEnergy');
  
end

function [y,energy,dynamic_sub] = RK4_FD_compact_spectral (u,deltat,N,mat_deriv1,mat_deriv2,F,h,constant_sub,factors_deriv,ind,filter)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with compact finite difference schemes
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
  
%%%%%% Get the Smagorinsky constant in case of dynamic model
  if (filter>0)
     kappa = 2; % filter ratio
     dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,factors_deriv,mat_deriv1);
     constant_sub = dynamic_sub ;
  else
     dynamic_sub = constant_sub ;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C = get_nonlinear_term(Un,ind(:,2:8),constant_sub,h,N,mat_deriv1,factors_deriv) ;
  
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un,ind(:,2:8),factors_deriv,mat_deriv2);
  
  k1  = C + vec_deriv2 ;
  Un2 = Un + deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C = get_nonlinear_term(Un2,ind(:,2:8),constant_sub,h,N,mat_deriv1,factors_deriv) ;
    
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un2,ind(:,2:8),factors_deriv,mat_deriv2);
  
  k2  = C + vec_deriv2 ;
  Un3 = Un + deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C = get_nonlinear_term(Un3,ind(:,2:8),constant_sub,h,N,mat_deriv1,factors_deriv) ;
    
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un3,ind(:,2:8),factors_deriv,mat_deriv2);

  k3  = C + vec_deriv2 ;
  Un4 = Un + deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C = get_nonlinear_term(Un4,ind(:,2:8),constant_sub,h,N,mat_deriv1,factors_deriv) ;
  
% Second derivative for the viscous term
  vec_deriv2 = get_second_derivative(Un4,ind(:,2:8),factors_deriv,mat_deriv2);
    
  k4 = C + vec_deriv2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  y = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 + F;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    energy = 0;
%    energy = get_energy(y,ind(:,2:8),mat_deriv1,mat_deriv2,factors_deriv,h,F);

endfunction

function vecC = get_nonlinear_term(Un,ind,constant_sub,h,N,mat_deriv1,factors_deriv)
% i-3  i-2  i-1  i  i+1  i+2  i+3 
%  1    2    3   4   5    6    7   
  
% Convective term, first the divergence form
   vec_derivUsquare = get_first_derivative( 0.5*Un.*Un , ind, factors_deriv,mat_deriv1 );
   
% Convective term, then the advective form
   vec_derivU = get_first_derivative( Un, ind, factors_deriv,mat_deriv1 );
   nonlin_adv = Un .* vec_derivU ;
   
% Convective term in skew-symmetric form
   vecC = - ( nonlin_adv + 2*vec_derivUsquare ) / 3 ;
   
% Subgrid term
   if (constant_sub>0)     
     vec_derivSG_im1 = vec_derivU(ind(:,3)) ; % derivative du/dx at node i-1
     vec_derivSG_im2 = vec_derivU(ind(:,2)) ; % derivative du/dx at node i-2
     vec_derivSG_im3 = vec_derivU(ind(:,1)) ; % derivative du/dx at node i-3
     vec_derivSG_ip1 = vec_derivU(ind(:,5)) ; % derivative du/dx at node i+1
     vec_derivSG_ip2 = vec_derivU(ind(:,6)) ; % derivative du/dx at node i+2
     vec_derivSG_ip3 = vec_derivU(ind(:,7)) ; % derivative du/dx at node i+3
     
     vec_derivSG = mat_deriv1 \ ( factors_deriv(5,1) * ( vec_derivSG_ip3.*abs(vec_derivSG_ip3) - vec_derivSG_im3.*abs(vec_derivSG_im3) ) + ...
                                  factors_deriv(4,1) * ( vec_derivSG_ip2.*abs(vec_derivSG_ip2) - vec_derivSG_im2.*abs(vec_derivSG_im2) ) + ...
                                  factors_deriv(3,1) * ( vec_derivSG_ip1.*abs(vec_derivSG_ip1) - vec_derivSG_im1.*abs(vec_derivSG_im1) ) );
     
     vecC += ( constant_sub * h )^2 * vec_derivSG ;
   endif
   
endfunction

function dfdx = get_first_derivative(Un,ind,factors_deriv,mat_deriv)
% i-3  i-2  i-1  i  i+1  i+2  i+3
%  1    2    3   4   5    6    7

% Compute the first spatial derivative with a compact finite difference scheme
  dfdx = factors_deriv(5,1) * ( Un(ind(:,7)) - Un(ind(:,1)) ) + ...
         factors_deriv(4,1) * ( Un(ind(:,6)) - Un(ind(:,2)) ) + ...
         factors_deriv(3,1) * ( Un(ind(:,5)) - Un(ind(:,3)) ) ;
  dfdx = mat_deriv \ dfdx ;
endfunction

function d2fdx2 = get_second_derivative(Un,ind,factors_deriv,mat_deriv)
% Compute the second spatial derivative with a compact finite difference scheme
  d2fdx2 = factors_deriv(5,2) * ( Un(ind(:,7)) -2*Un(ind(:,4)) + Un(ind(:,1)) ) + ...
           factors_deriv(4,2) * ( Un(ind(:,6)) -2*Un(ind(:,4)) + Un(ind(:,2)) ) + ...
           factors_deriv(3,2) * ( Un(ind(:,5)) -2*Un(ind(:,4)) + Un(ind(:,3)) ) ; 
  d2fdx2 = mat_deriv \ d2fdx2;
endfunction

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

function dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter,factors_deriv,mat_deriv)
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
##   u_filter   = apply_filter(Un    ,ind,filter) ;
##   usq_filter = apply_filter(Un.*Un,ind,filter) ;
##   u_filtersq = u_filter.*u_filter;
##   
##   H = get_first_derivative( u_filtersq*0.5 , ind(:,2:8) ,factors_deriv,mat_deriv) - ...
##       get_first_derivative( usq_filter*0.5 , ind(:,2:8) ,factors_deriv,mat_deriv ) ;
##   
##   deriv_u_filter = get_first_derivative(u_filter,ind(:,2:8),factors_deriv,mat_deriv);
##   tmp1 = apply_filter(deriv_u_filter.*abs(deriv_u_filter),ind,filter);
##   
##   M = kappa*kappa*get_first_derivative( abs(deriv_u_filter).*deriv_u_filter , ind(:,2:8) ,factors_deriv,mat_deriv ) -...
##       get_first_derivative( tmp1 , ind(:,2:8) ,factors_deriv,mat_deriv );
##   
##   csdsq = sum(H.*M) / sum(M.*M); % (Cs * Delta)^2
##   dynamic_sub = sqrt(abs(csdsq)) / h ;
   
% New method
   u_filter = apply_filter(Un ,ind,filter) ;
   L        = apply_filter(Un.*Un ,ind,filter) - u_filter.*u_filter ;
   
   deriv_u  =  get_first_derivative(Un,ind(:,2:8),factors_deriv,mat_deriv);
   deriv_u_filter = apply_filter(deriv_u,ind,filter);
   M = kappa*kappa* deriv_u_filter.*abs(deriv_u_filter) - ...
       apply_filter( deriv_u .*abs(deriv_u) ,ind,filter) ;
   
   csdsq = 0.5 * sum(L.*M) / sum(M.*M); % (Cs * Delta)^2
   dynamic_sub = sqrt(abs(csdsq)) / h ;
endfunction

function energy = get_energy(Un,ind,mat_deriv1,mat_deriv2,factors_deriv,h,F)
% Compute the numerical energy produced by the spatial discretization of the 
% convective term in the skew-symmetric form
  energy = 0;
  
% Divergence form
  deriv1_u = get_first_derivative( 0.5*Un.*Un ,ind,factors_deriv,mat_deriv1);
  energy1 = h * Un' * deriv1_u;
    
% Advective form
  deriv1_u = get_first_derivative(Un,ind,factors_deriv,mat_deriv1);
  energy2 = h * (Un .* Un)' * deriv1_u; 

% Sum between advective and divergent forms for nonlinear term
  energy = energy - (2*energy1 + energy2)/3 ;
  
% % Contribution for viscous dissipation
%  vec_deriv2 = get_second_derivative(Un,ind,factors_deriv,mat_deriv2) ;
%  energy = energy + h * Un' * deriv2_u;
     
% % Contribution of forcing term
%  energy = energy + h * Un' * F;
endfunction