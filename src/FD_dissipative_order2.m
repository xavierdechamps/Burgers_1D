function FD_dissipative_order2 (N,nu,constant_sub,filter,L,time,nbrpointtemp,name,file_spectrum)
% Solve the 1D forced Burgers equation with the energy dissipative Hd2 scheme of 
% order 2 for the convective term
% The unknown of the equation is the velocity, thus 1 unknown per node
%
% du     du     d2u
% --  + u-- = nu--- + F
% dt     dx     dx2
%************* Initialization of the parameters **************************
  disp("*********************************************************************")  
  disp("Finite difference - 2nd order divergence form for convective term ")
  disp("*********************************************************************")
  
  h=L/(N);% Length of the elements
  X = linspace(0,L,N)';
  nbrpointtime = nbrpointtemp;
  deltat       = time/nbrpointtime;

% ************* Initial condition on the solution ************************
% Random solution for turbulent flow
  u(:,1)=2*rand(N,1)-1;
% Sinus solution for non-forced Burgers equation
%  u(:,1) = sin(X * 2*pi/L);
%  u(:,1)=0.5 + sin(X(:)* 2*pi/L);

  timeBeforeStatistics = 10;

  kinEnergy = zeros(nbrpointtemp+1,1); kinEnergy(1) = u(:,1)' * u(:,1) * h * 0.5;
%  filename=[name,num2str(1),'.mat']; uu=u(:,1); save(filename,'uu');

  z=2;j=2;
  nbrPointsStatistics=0; kinEnergyMean=0;
  
  spectralEnergy=zeros(N,1);
  reference_spectrum=load(file_spectrum);
  
% Store the indices of the neighbour nodes i-2 i-1 i i+1 i+2 used in the derivatives
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
  diverged = false;
  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand(); phi3=2*pi*rand();
    KForce=2*pi/L;
    F=sqrt(deltat)*(cos(2*KForce*X+phi2)+cos(3*KForce*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F=0;

%******** Call Runge-Kutta and compute kinematic energy ********
    [u(:,z),dynamic_smag_constant(i-1)] = ...
              RK4_FD_dissipative_order2 (u(:,z-1),deltat,N,nu,F,h,constant_sub,ind,filter);

    kinEnergy(i) = h*0.5*u(:,z)'*u(:,z);
    
    if (i*deltat>=timeBeforeStatistics)
%      kinEnergyMean = kinEnergyMean*nbrPointsStatistics/(nbrPointsStatistics+1) + kinEnergy(i)/(nbrPointsStatistics+1);  
      fft_u = fft(u(:,end));
      newSpectralEnergy = fft_u.*conj(fft_u)/(N)/(N);
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
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done - Statistics are stored' ))
%           filename2=strcat('Spectral_energy_',name,'.mat'); spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics; save(filename2,'spectralEnergyOut');
        else
           disp(strcat(num2str(100*i/(nbrpointtime+1) ),' % done' ))
        end
       
        subplot(2,2,1)
        plot([X;X(end)+h]/L, [uu; uu(1)],'Linewidth',3)
        grid on; xlabel('x/(2*\pi)'); ylabel('u(t)')
        xlim([0 1])
        title(strcat('Time= ',num2str(i*deltat),', Re= ',num2str(mean(uu)*L/nu)))
        
        subplot(2,2,3)
        plot((0:(i-1))*deltat,kinEnergy(1:i),'Linewidth',3)
        grid on; xlabel('Time'); ylabel('E(t)')
        xlim([0 time])
         
        subplot(2,2,2)
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
        grid on; xlabel('k'); ylabel('E(k)')
        xlim([1 reference_spectrum(end,1)])
         
        subplot(2,2,4)
        plot((1:(i-1))*deltat,dynamic_smag_constant(1:i-1),'Linewidth',3)
        grid on; xlabel('Time'); ylabel('Smagorinsky C_s(t)') ;
        xlim([0 time])
        
        drawnow
    end
    
    z=z+1;
  
    CFL=u(:,end)*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(['Divergence of ',name,', CFL=',num2str(max(CFL))]);
        diverged = true;
        break;
    end
  end
  
  mean_Smagorinsky = mean(dynamic_smag_constant)

  spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  if (~diverged)
    save(filename2,'-ascii','spectralEnergyOut');
  end
  
  %filename=strcat('Energy_',name,'.mat');
  %save(filename,'kinEnergy');

end

function [y,dynamic_sub] = RK4_FD_dissipative_order2 (u,deltat,N,nu,F,h,constant_sub,ind,filter)
% Temporal integration of the 1D Burgers equation with an explicit 4 steps Runge-Kutta scheme
% Spatial discretization with an energy dissipative Hd2 scheme of 
% order 2 for the convective term
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
  
  if (filter>0)
     kappa = 2; % filter ratio
     dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter);
     constant_sub = dynamic_sub ;
  else
     dynamic_sub = constant_sub ;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C      = get_nonlinear_term(Un,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un , ind(:,4:6) , h ) ;
	k1     = du2dx2 + C;
	Un2    = Un+deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C      = get_nonlinear_term(Un2,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un2 , ind(:,4:6) , h ) ;
	k2     = du2dx2 + C;
	Un3    = Un+deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
	C      = get_nonlinear_term(Un3,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un3 , ind(:,4:6) , h ) ;
	k3     = du2dx2 + C;
	Un4    = Un+deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
	C      = get_nonlinear_term(Un4,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un4 , ind(:,4:6) , h ) ;
	k4     = du2dx2 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y = Un + deltat*(k1 + 2*k2 + 2*k3 +k4 )/6 +F;
  
end

function vecC = get_nonlinear_term(Un,ind,constant_sub,h,N)
% Compute the non-linear + subgrid terms with the dissipative scheme
% i-4  i-3  i-2  i-1  i  i+1  i+2  i+3  i+4
%  1    2    3    4   5   6    7    8    9
  
% non-linear term in divergence form
   du2dx = get_first_derivative( Un .* Un , ind(:,4:6) , h );
   vecC  = - du2dx * 0.5 ; 
   
% Subgrid term
   if (constant_sub>0)
     dudx_im1 = get_first_derivative( Un       , ind(:,3:5) , h ); % derivative at node i-1
     dudx_ip1 = get_first_derivative( Un       , ind(:,5:7) , h ); % derivative at node i+1
     vecC    += constant_sub^2 * h * ( dudx_ip1.*abs(dudx_ip1) - dudx_im1.*abs(dudx_im1) ) * 0.5;
   endif

endfunction

function dudx = get_first_derivative(Un,ind,h)
% Un has size (N,3)
   dudx = ( Un(ind(:,3)) - Un(ind(:,1)) ) * 0.5 / h ;
endfunction

function du2dx2 = get_second_derivative(Un,ind,h)
% Un has size (N,3)
   du2dx2 = ( Un(ind(:,1)) - 2*Un(ind(:,2)) + Un(ind(:,3)) ) / (h*h) ;
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

function dynamic_sub = get_dynamic_smagorinsky(Un,ind,h,kappa,filter)
% Compute the Smagorinsky constant by a dynamic model
% See "Evaluation of explicit and implicit LES closures for Burgers turbulence"
% by R. Maulik and O. San, Journal of Computational and Applied Mathematics 327 (2018) 12-40
   u_filter   = apply_filter(Un    ,ind,filter) ;
   usq_filter = apply_filter(Un.*Un,ind,filter) ;
   u_filtersq = u_filter.*u_filter;
   
   H = get_first_derivative( u_filtersq*0.5 , ind(:,4:6) , h) - ...
       get_first_derivative( usq_filter*0.5 , ind(:,4:6) , h ) ;
   
   deriv_u_filter = get_first_derivative(u_filter,ind(:,4:6),h);
   tmp1 = apply_filter(deriv_u_filter.*abs(deriv_u_filter),ind,filter);
   
   M = kappa*kappa*get_first_derivative( abs(deriv_u_filter).*deriv_u_filter , ind(:,4:6) , h ) -...
       get_first_derivative( tmp1 , ind(:,4:6) , h );
   
   csdsq = sum(H.*M) / sum(M.*M); % (Cs * Delta)^2
   dynamic_sub = sqrt(abs(csdsq)) / h ;
endfunction