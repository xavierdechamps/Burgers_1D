function energy = get_kinematic_energy(h,NG,u,N,type)
% Integrate ( u^2 ) / 2 over the whole domain with numerical integration
w=zeros(NG); ksi=zeros(NG);
fourth = 0.25;
half   = 0.5;
eigth  = 0.125;

if (NG==1) 
    w(1)=2; 
    ksi(1)=0;
elseif (NG==2) 
    w(:)=1; 
    ksi(1)=0.577350269189626; ksi(2)=-ksi(1);
elseif (NG==3) 
    w(1)=0.555555555555556; w(2)=w(1); w(3)=0.888888888888889;
    ksi(1)=0.774596669241483; ksi(2)=-ksi(1); ksi(3)=0;
elseif(NG==4) 
    w(1)=0.347854845137454; w(2)=w(1); w(3)=0.652145154862546; w(4)=w(3);
    ksi(1)=0.861136311594953; ksi(2)=-ksi(1); ksi(3)=0.339981043584856; ksi(4)=-ksi(3);
end

energy=0;

if (type==1) % Lagrange linear elements

   ind = zeros(N,2);
   ind(:,1)   = 1:N; 
   ind(:,2)   = 2:N+1; 
   ind(N,:) = [N 1];
   
   for j=1:NG % Sommer sur tous les points de Gauss
       energy = energy + w(j)* ( -(ksi(j)-1) * u(ind(:,1))' + (ksi(j)+1) * u(ind(:,2))' ) ...
                             * ( -(ksi(j)-1) * u(ind(:,1))  + (ksi(j)+1) * u(ind(:,2))  );
   end
   energy = energy * h * eigth * half;
    
elseif (type==2) % Cubic Hermite elements
    for i=1:2:2*N-3 % sum over all the elements
        if (i==2*N-3) 
            a=[2*N-3 1];
        else
            a=[i i+2];
        end
        for j=1:NG % sum over all the Gauss points
            energy=energy + h*fourth*w(j) * (u(a(1))*(1+fourth*(ksi(j)-2)*(ksi(j)+1)^2) + u(a(1)+1)*h*fourth*(ksi(j)+1)*(-ksi(j)+fourth*(ksi(j)+1)^2) +...
                u(a(2))*fourth*(2-ksi(j))*(ksi(j)+1)^2 + u(a(2)+1)*half*h*eigth*(ksi(j)-1)*(ksi(j)+1)^2)^2;
        end
    end
end
