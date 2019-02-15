function FD_nonlinear_schemes (N,nu,constant_sub,L,time,nbrpointtemp,name,file_spectrum,submethod,order_visc)
% Solve the 1D forced Burgers equation with a nonlinear discretization for the convective term
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("*********************************************************************")  
  disp("Finite difference - nonlinear scheme for convective term")
  disp("*********************************************************************")  

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

  kinEnergy    = zeros(nbrpointtemp+1,1); kinEnergy(1) = u(:,1)' * u(:,1) * h * 0.5;
  
%  filename=[name,num2str(1),'.mat']; uu=u(:,1); save(filename,'uu');

  z=2;j=2; ind_error=1;
  nbrPointsStatistics=0; kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);
  
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
  
  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand();    phi3=2*pi*rand();
    KForce=2*pi/L;
    F = sqrt(deltat)*(cos(2*KForce*X+phi2)+cos(3*KForce*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F = 0;
    
%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z) = RK4_FD_nonlinear(u(:,z-1),deltat,N,nu,F,h,constant_sub,ind,submethod,order_visc);
    
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
%       filename=[name,num2str(j),'.mat']; save(filename,'uu'); j=j+1;
        clear u;
        z=1;
        u(:,1)=uu;
        
        reynolds_number        = mean(uu)*L/nu;
        kolmogorov_length      = L*reynolds_number^(-0.75);
        kolmogorov_wave_length = 2*pi/kolmogorov_length ;
                
        if (nbrPointsStatistics > 0)
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
        
        subplot(2,2,1)
        plot([X; X(end)+h]/L, [uu; uu(1)],'b','Linewidth',3)
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
         
         subplot(1,2,2)
         loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
         grid on; xlabel('k'); ylabel('E(k)')
         xlim([1 reference_spectrum(end,1)])
        
         drawnow
    end
    
    z=z+1;
  
    CFL=u(:,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(['Divergence of ',name]);
        break;
    end
  end
  
%  relative_error
  
  spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  save(filename2,'spectralEnergyOut');
  
  %filename=strcat('Energy_',name,'.mat');
  %save(filename,'kinEnergy');
  
end

function y = RK4_FD_nonlinear(u,deltat,N,nu,F,h,constant_sub,ind,submethod,order_visc)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with nonlinear discretizations for the convective term
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
  Un=u;

  switch submethod
       case 1
           function_nonlinear = 'get_nonlinear_term_UP1' ;
       case 2
           function_nonlinear = 'get_nonlinear_term_UP2' ;
       case 3
           function_nonlinear = 'get_nonlinear_term_UP3' ;
       case 4
           function_nonlinear = 'get_nonlinear_term_DFD' ;
       case 5
           function_nonlinear = 'get_nonlinear_term_TVD' ;
       otherwise
           disp('Unknown kind of nonlinear scheme, exiting the code...')
           return
  end
%% WENO IS NOT FULLY FUNCTIONING YET
%    u(:,z) = RK4_FD_WENO5 (u(:,z-1),deltat,N,K,F,h,constant_sub);
%    u(:,z) = NSSP_RK5_FD_WENO5 (u(:,z-1),deltat,N,K,F,h,constant_sub);
  
  switch order_visc
       case 2
           function_visc = 'get_second_derivative_order2' ;
       case 4
           function_visc = 'get_second_derivative_order4' ;
       case 6
           function_visc = 'get_second_derivative_order6' ;
       otherwise
           disp('Unknown kind of viscous scheme, exiting the code...')
           return
  end

  function_handle_nonlinear = str2func(function_nonlinear) ;
  function_handle_viscous   = str2func(function_visc) ;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C      = function_handle_nonlinear(Un,ind,h) ;
  du2dx2 = nu * function_handle_viscous( Un , ind , h ) ;
  k1     = du2dx2 + C;
  Un2    = Un + deltat * 0.5 * k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C      = function_handle_nonlinear(Un2,ind,h) ;
  du2dx2 = nu * function_handle_viscous( Un2 , ind , h ) ;
  k2     = du2dx2 + C;
  Un3    = Un + deltat * 0.5 * k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
  C      = function_handle_nonlinear(Un3,ind,h) ;
  du2dx2 = nu * function_handle_viscous( Un3 , ind , h ) ;
  k3     = du2dx2 + C;
  Un4    = Un + deltat * k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
  C      = function_handle_nonlinear(Un4,ind,h) ;
  du2dx2 = nu * function_handle_viscous( Un4 , ind , h ) ;
  k4     = du2dx2 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y  = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 +F;
end

function vecC = get_nonlinear_term_UP1(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9

  vpi = max( Un(ind(:,5)) ,0); % v+ at node i
  vmi = min( Un(ind(:,5)) ,0); % v- at node i
 
  vpip1 = max( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v+ at node i+1/2
  vpim1 = max( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v+ at node i-1/2
  vmip1 = min( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v- at node i+1/2
  vmim1 = min( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v- at node i-1/2
  
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = - ( ( vpi + vpip1 - vmi - vmim1 ) .* Un(ind(:,5)) + ...
             ( vmi + vmip1 ) .* Un(ind(:,6)) - ...
             ( vpi + vpim1 ) .* Un(ind(:,4)) ) / ( 3. * h ) ;
endfunction

function vecC = get_nonlinear_term_UP2(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9

  vpi = max( Un(ind(:,5),1) ,0); % v+ at node i
  vmi = min( Un(ind(:,5),1) ,0); % v- at node i
  
  vpip1 = max( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v+ at node i+1/2
  vpim1 = max( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v+ at node i-1/2
  vmip1 = min( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v- at node i+1/2
  vmim1 = min( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v- at node i-1/2
  
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = - ( ( vpi + vpip1 ) .* ( 3*Un(ind(:,5)) - Un(ind(:,4)) ) + ...
             ( vmi + vmip1 ) .* ( 3*Un(ind(:,6)) - Un(ind(:,7)) ) - ...
             ( vpi + vpim1 ) .* ( 3*Un(ind(:,4)) - Un(ind(:,3)) ) - ...
             ( vmi + vmim1 ) .* ( 3*Un(ind(:,5)) - Un(ind(:,6)) ) ) / ( 6 * h ) ;
endfunction

function vecC = get_nonlinear_term_UP3(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9

  vpi = max( Un(ind(:,5),1) ,0); % v+ at node i
  vmi = min( Un(ind(:,5),1) ,0); % v- at node i
  
  vpip1 = max( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v+ at node i+1/2
  vpim1 = max( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v+ at node i-1/2
  vmip1 = min( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v- at node i+1/2
  vmim1 = min( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v- at node i-1/2
  
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = - ( ( vpi + vpip1 ) .* ( 2*Un(ind(:,6)) + 5*Un(ind(:,5)) - Un(ind(:,4)) ) + ...
             ( vmi + vmip1 ) .* ( 2*Un(ind(:,5)) + 5*Un(ind(:,6)) - Un(ind(:,7)) ) - ...
             ( vpi + vpim1 ) .* ( 2*Un(ind(:,5)) + 5*Un(ind(:,4)) - Un(ind(:,3)) ) - ...
             ( vmi + vmim1 ) .* ( 2*Un(ind(:,4)) + 5*Un(ind(:,5)) - Un(ind(:,6)) ) ) / (18*h);
endfunction

function vecC = get_nonlinear_term_DFD(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
  f = 0.21 ;
  
  indip1 = ind(:,4:8);
  indim1 = ind(:,2:6);
  
  cip1 = - 1 ./ ( 6 * (1 + f * max( min( ( Un(indip1(:,5)) - 4*Un(indip1(:,4)) + 6*Un(indip1(:,3)) - 4*Un(indip1(:,2)) + Un(indip1(:,1)) ) ./ ( Un(indip1(:,4)) - 2*Un(indip1(:,3)) + Un(indip1(:,2)) ) , 0 ) , -3 ) ) );
  cim1 = - 1 ./ ( 6 * (1 + f * max( min( ( Un(indim1(:,5)) - 4*Un(indim1(:,4)) + 6*Un(indim1(:,3)) - 4*Un(indim1(:,2)) + Un(indim1(:,1)) ) ./ ( Un(indim1(:,4)) - 2*Un(indim1(:,3)) + Un(indim1(:,2)) ) , 0 ) , -3 ) ) );
    
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = - ( Un(ind(:,5)) .* ( Un(ind(:,6)) + cip1.*( Un(ind(:,7)) - 2*Un(ind(:,6)) + Un(ind(:,5)) ) - ...
                               Un(ind(:,4)) - cim1.*( Un(ind(:,5)) - 2*Un(ind(:,4)) + Un(ind(:,3)) ) ) + ... 
             Un(ind(:,6)).^2 + cip1 .* ( Un(ind(:,7)).^2 -2*Un(ind(:,6)).^2 + Un(ind(:,5)).^2 ) - ...
             Un(ind(:,4)).^2 - cim1 .* ( Un(ind(:,5)).^2 -2*Un(ind(:,4)).^2 + Un(ind(:,3)).^2 ) ) / ( 6 * h ) ;
endfunction

function vecC = get_nonlinear_term_TVD(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9

  vpi = max( Un(ind(:,5),1) ,0); % v+ at node i
  vmi = min( Un(ind(:,5),1) ,0); % v- at node i
  
  vpip1 = max( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v+ at node i+1/2
  vpim1 = max( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v+ at node i-1/2
  vmip1 = min( ( Un(ind(:,5)) + Un(ind(:,6)) ) * 0.5   ,0); % v- at node i+1/2
  vmim1 = min( ( Un(ind(:,5)) + Un(ind(:,4)) ) * 0.5   ,0); % v- at node i-1/2
  
  rim1 = ( Un(ind(:,4),1) - Un(ind(:,3),1) ) ./ ( Un(ind(:,5),1) - Un(ind(:,4),1) );
  ri   = ( Un(ind(:,5),1) - Un(ind(:,4),1) ) ./ ( Un(ind(:,6),1) - Un(ind(:,5),1) );
  rip1 = ( Un(ind(:,6),1) - Un(ind(:,5),1) ) ./ ( Un(ind(:,7),1) - Un(ind(:,6),1) );
  
% Skew-symmetric form = energy conservative form (see the appendix in Fauconnier - Dick)
  vecC = - ( ( vpi + vpip1 ) .* ( Un(ind(:,5),1) + 0.5*get_psi_TVD(ri)     .*(Un(ind(:,6),1)-Un(ind(:,5),1)) ) - ...
             ( vpi + vpim1 ) .* ( Un(ind(:,4),1) + 0.5*get_psi_TVD(rim1)   .*(Un(ind(:,5),1)-Un(ind(:,4),1)) ) + ...
             ( vmi + vmip1 ) .* ( Un(ind(:,6),1) - 0.5*get_psi_TVD(1./rip1).*(Un(ind(:,6),1)-Un(ind(:,5),1)) ) - ...
             ( vmi + vmim1 ) .* ( Un(ind(:,5),1) - 0.5*get_psi_TVD(1./ri)  .*(Un(ind(:,5),1)-Un(ind(:,4),1)) ) ) / ( 3. * h ) ;
endfunction

function psi = get_psi_TVD(r)
  psi = max ( min ( r , 2 ) , 0 );
endfunction

function du2dx2 = get_second_derivative_order2(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
   du2dx2 = ( Un(ind(:,4)) - 2*Un(ind(:,5)) + Un(ind(:,6)) ) / (h*h) ;
endfunction

function du2dx2 = get_second_derivative_order4(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
   du2dx2 = ( - Un(ind(:,3)) + 16*Un(ind(:,4)) - 30*Un(ind(:,5)) + ...
                               16*Un(ind(:,6)) - Un(ind(:,7)) ) / (12*h*h) ;
endfunction

function du2dx2 = get_second_derivative_order6(Un,ind,h)
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
   du2dx2 = ( Un(ind(:,2))/90 - 3*Un(ind(:,3))/20 + 1.5*Un(ind(:,4)) - ...
              49*Un(ind(:,5))/18 + 1.5*Un(ind(:,6)) - 3*Un(ind(:,7))/20 + ...
              Un(ind(:,8))/90 ) / (h*h) ;
endfunction