function [ one,two,three,four,five ] = Calc_Coefficient( n,DR,DQ,dt,dQ,dr,Q,r,one,two,three,four,five )
% function that calc jacobi coefficient for the heat equation
% recive grid information
one(1:n,1:n) = 3+(4*dt)./(DR(1:n,1:n).*dr)+(4*dt)./((DQ(1:n,1:n).*dQ).*r(1:n,1:n).^2);
two(1:n,1:n) = (4*dt)./(r(1:n,1:n).*(DR(1:n,1:n)+dr)) +(2*dt)./(DR(1:n,1:n).*dr);
three(1:n,1:n) = -(4*dt)./(r(1:n,1:n).*(DR(1:n,1:n)+dr)) +(2*dt)./(DR(1:n,1:n).*dr);
four(1:n,1:n) = (2*dt)./((r(1:n,1:n).^2).*(tan(Q(1:n,1:n)).*(DQ(1:n,1:n)+dQ))) + (2*dt)./((DQ(1:n,1:n).*dQ).*r(1:n,1:n).^2);
five(1:n,1:n) = -(2*dt)./((r(1:n,1:n).^2).*(tan(Q(1:n,1:n)).*(DQ(1:n,1:n)+dQ))) + (2*dt)./((DQ(1:n,1:n).*dQ).*r(1:n,1:n).^2);

end

