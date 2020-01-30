function [ RHS ] = CreateRHS( n ,Tn,F_T)
% this function creates RHS vector for solving Newton-Raphson Method
% Recive n,Tn,F_T 
RHS(1:n^2)=NaN;

RHS(n+1:n^2-2)=F_T(n+1:n^2-2);
RHS(2:n-1)=(-3*Tn(2:n-1)-1*Tn(2+2*n:3*n-1)+4*Tn(2+n:2*n-1));  % n-2 start dt/dq=0
RHS(n^2-n+2:n^2-1)=(-3*Tn(n^2-n+2:n^2-1)-1*Tn(n^2-3*n+2:n^2-2*n-1)+4*Tn(n^2-2*n+2:n^2-n-1)); % n-2 end dt/dq=0
RHS(1:n:n^2)=Tn(1:n:n^2);
RHS(n:n:n^2)=Tn(n:n:n^2)-1;
RHS=RHS';

end

