function sol = get_analytical_solution(t,nu,x,n)

	A0 = besseli(0,1/(2*pi*nu)) ;

	denom = A0 ;
	num   = 0;

	for k = 1:n
		Ai = besseli(k,1/(2*pi*nu)) ;

		addNum   = k * Ai * exp(-k*k*pi*pi*t*nu) * sin(k*pi*x) ;
		addDenom = 2 * Ai * exp(-k*k*pi*pi*t*nu) * cos(k*pi*x) ;

		num   = num   + addNum;
		denom = denom + addDenom;

	end

	sol = 4*pi*nu*num./denom ;

end
