function y = Interpolation(u,h,N,polyn,nbr_points)
% High order finite element schemes show oscillations within the elements. 
% To put forward this phenomenon, an interpolation of the solution is made within
% the elements with nbr_points points.
		
	y = zeros(nbr_points*N,1);
		
    if (polyn == 1) % cubic Hermite
    
		ind = zeros(N,4);
		ind(:,1) = 1:2:2*(N-1)+1;
		ind(:,2) = 2:2:2*(N-1)+2;
		ind(:,3) = 3:2:2*(N-1)+3;
		ind(:,4) = 4:2:2*(N-1)+4;
		ind(N,:) = [2*N-1 2*N 1 2];

		fct_in_element = eval_hermite_p3(nbr_points,h) ;
		
		for i=1:nbr_points
			y(i:nbr_points:end) = u(ind(:,1))*fct_in_element(1,i) + ...
                            u(ind(:,2))*fct_in_element(2,i) + ...
                            u(ind(:,3))*fct_in_element(3,i) + ...
                            u(ind(:,4))*fct_in_element(4,i);
		end
		
	elseif (polyn == 2) % cubic Lagrange P3
		ind = zeros(N,4);
		ind(:,1) = 1:3:3*N-2;
		ind(:,2) = 2:3:3*N-1;
		ind(:,3) = 3:3:3*N;
		ind(:,4) = 4:3:3*N+1;
		ind(N,:) = [3*N-2 3*N-1 3*N 1];

		fct_in_element = eval_lagrange_p3(nbr_points);
		
		for i=1:nbr_points
			y(i:nbr_points:end) = u(ind(:,1))*fct_in_element(1,i) + ...
                            u(ind(:,2))*fct_in_element(2,i) + ...
                            u(ind(:,3))*fct_in_element(3,i) + ...
                            u(ind(:,4))*fct_in_element(4,i);
		end 
		
	elseif (polyn == 3) % pentatonic Hermite
		ind = zeros(N,6);
		ind(:,1) = 1:3:3*N-2; 
		ind(:,2) = 2:3:3*N-1; 
		ind(:,3) = 3:3:3*N;
		ind(:,4) = 4:3:3*N+1;
		ind(:,5) = 5:3:3*N+2;
		ind(:,6) = 6:3:3*N+3;
		ind(end,4:6) = 1:3;

		fct_in_element = eval_hermite_p5(nbr_points,h) ;
		
		for i=1:nbr_points
			y(i:nbr_points:end) = u(ind(:,1))*fct_in_element(1,i) + ...
                            u(ind(:,2))*fct_in_element(2,i) + ...
                            u(ind(:,3))*fct_in_element(3,i) + ...
                            u(ind(:,4))*fct_in_element(4,i) + ...
                            u(ind(:,5))*fct_in_element(5,i) + ...
                            u(ind(:,6))*fct_in_element(6,i);
		end
		
	else
		error('Unknown kind of interpolation')
	end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fct_in_element = eval_hermite_p3(n,h)
	x = ( 0:1:(n-1) )/n;
	fct_in_element = zeros(4,n);
		
	fct_in_element(1:4,1:n) = [ 1 - 3*x.*x + 2*x.*x.*x       ; 
						                  h * ( x - 2*x.*x + x.*x.*x ) ; 
						                  3*x.*x - 2*x.*x.*x           ;  
						                  h * (x.*x.*x - x.*x)
						                ] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fct_in_element = eval_lagrange_p3(n)
	x = ( 0:1:(n-1) )/n;
	fct_in_element = zeros(4,n);
	
	fct_in_element(1:4,1:n) = [ (1-x).*(2-3*x).*(1-3*x)*0.5  ; 
						                  4.5*x.*(1-x).*(2-3*x)        ; 
						                  4.5*x.*(1-x).*(3*x-1)        ;  
						                  x.*(2-3*x).*(1-3*x)*0.5 
						                ] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fct_in_element = eval_hermite_p5(n,h)
	x = ( 0:1:(n-1) )/n;
	fct_in_element = zeros(6,n);
			
	fct_in_element(1:6,1:n) = [ 1 - 10*x.^3 + 15*x.^4 - 6*x.^5          ; 
						                  h * ( x -6*x.^3 + 8*x.^4 - 3*x.^5  )    ; 
						                  h*h*(x.^2 - 3*x.^3 + 3*x.^4 - x.^5)*0.5 ;  
						                  10*x.^3 - 15*x.^4 + 6*x.^5              ;
						                  h*(-4*x.^3 + 7*x.^4 - 3*x.^5 )          ;
						                  h*h* (x.^3 - 2*x.^4 + x.^5)*0.5
						                ] ;
end
