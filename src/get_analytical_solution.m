function sol = get_analytical_solution(t,nu,x,n)

	A0 = besseli(0,1/(2*pi*nu)) ;
	
	denom = A0 ;
	num   = 0;

	for i = 1:n	
		Ai = besseli(i,1/(2*pi*nu)) ;

		addNum   = i * Ai * exp(-i*i*pi*pi*t*nu) * sin(i*pi*x) ;
		addDenom = 2 * Ai * exp(-i*i*pi*pi*t*nu) * cos(i*pi*x) ;
		
		num   = num   + addNum;
		denom = denom + addDenom;
	
	end

	sol = 4*pi*nu*num./denom ;

end
