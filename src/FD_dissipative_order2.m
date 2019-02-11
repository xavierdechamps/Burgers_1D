function FD_dissipative_order2 (N,nu,constant_sub,L,time,nbrpointtemp,name,file_spectrum)
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
  
  for i=2:nbrpointtime+1   
%***************** Forcing term with with-noise random phase ******************
    phi2=2*pi*rand(); phi3=2*pi*rand();
    KForce=2*pi/L;
    F=sqrt(deltat)*(cos(2*KForce*X+phi2)+cos(3*KForce*X+phi3));
% Uncomment the following line in case of sinus wave with non-linear convection
%    F=0;

%******** Call Runge-Kutta and compute kinematic energy ********
    u(:,z) = RK4_FD_dissipative_order2 (u(:,z-1),deltat,N,nu,F,h,constant_sub);

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
         
        subplot(1,2,2)
        loglog(0:(N/2-1),spectralEnergy(1:(N/2))/nbrPointsStatistics,'r','Linewidth',3, reference_spectrum(:,1),reference_spectrum(:,2),'b','Linewidth',3)
        grid on; xlabel('k'); ylabel('E(k)')
         xlim([1 reference_spectrum(end,1)])
        
        drawnow
    end
    
    z=z+1;
  
    CFL=u(:,end).*deltat/h;
% Stability criterion for explicit Runge Kutta 4
    if (max(CFL)>2.8)
        disp(strcat('Divergence of ',name,', CFL=',num2str(max(CFL))));
        break;
    end
  end

  spectralEnergyOut = spectralEnergy(1:(N/2))/nbrPointsStatistics;
  filename2=strcat('Spectral_energy_',name,'.mat');
  save(filename2,'spectralEnergyOut');
  
  %filename=strcat('Energy_',name,'.mat');
  %save(filename,'kinEnergy');

end

function y= RK4_FD_dissipative_order2 (u,deltat,N,nu,F,h,constant_sub)
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

  ind = zeros(N,5);
  ind(:,3) = 1:N;
  ind(:,2) = circshift(ind(:,3), 1,1);
  ind(:,1) = circshift(ind(:,3), 2,1);
  ind(:,4) = circshift(ind(:,3),-1,1);
  ind(:,5) = circshift(ind(:,3),-2,1);
  
	sixth = 1/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second step
  C      = get_nonlinear_term(Un,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un(ind(:,2:4)) , h ) ;
	k1     = du2dx2 + C;
	Un2    = Un+deltat*0.5*k1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Third step
  C      = get_nonlinear_term(Un2,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un2(ind(:,2:4)) , h ) ;
	k2     = du2dx2 + C;
	Un3    = Un+deltat*0.5*k2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourth step
	C      = get_nonlinear_term(Un3,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un3(ind(:,2:4)) , h ) ;
	k3     = du2dx2 + C;
	Un4    = Un+deltat*k3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fifth step
	C      = get_nonlinear_term(Un4,ind,constant_sub,h,N);
  du2dx2 = nu * get_second_derivative( Un4(ind(:,2:4)) , h ) ;
	k4     = du2dx2 + C;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y = Un + deltat*sixth*(k1 + 2*k2 + 2*k3 +k4 ) +F;
  
end

function vecC = get_nonlinear_term(Un,ind,constant_sub,h,N)
% Compute the non-linear + subgrid terms with the dissipative scheme
% i = i-2  i-1  i  i+1  i+2
%      1    2   3   4    5
##	fourth = 0.25;
##	one_over_h = 1/h;
  
% non-linear terms
% Divergence form
   du2dx = get_first_derivative( Un(ind(:,2:4)) .* Un(ind(:,2:4)) , h );
   vecC  = - du2dx * 0.5 ; 
   
% Subgrid term
   if (constant_sub>0)
     dudx_im1 = get_first_derivative( Un(ind(:,1:3)) , h ); % derivative at node i-1
     dudx_ip1 = get_first_derivative( Un(ind(:,3:5)) , h ); % derivative at node i+1
     vecC    += constant_sub^2 * h * ( dudx_ip1.*abs(dudx_ip1) - dudx_im1.*abs(dudx_im1) ) * 0.5;
   endif

endfunction

function dudx = get_first_derivative(Un,h)
% Un has size (N,3)
   dudx = ( Un(:,3) - Un(:,1) ) * 0.5 / h ;
endfunction

function du2dx2 = get_second_derivative(Un,h)
% Un has size (N,3)
   du2dx2 = ( Un(:,1) - 2*Un(:,2) + Un(:,3) ) / (h*h) ;
endfunction
