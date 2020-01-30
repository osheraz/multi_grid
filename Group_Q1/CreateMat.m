function [ A ] = CreateMat( A,n ,DF_T_ij,DF_T_ijminus1,DF_T_ijplus1,DF_T_iminus1j,DF_T_iplus1j)
% this function get n (size for grid n^2*n^2) and the parameters for all
% diag values and return the Jacobian matrix A. 


%% center diag 
% A(a:b:c)=>a is first cell, b skipping between cells, c last cell
A((n+1)*n^2+n+2:n^2+1:n^4-n^3-2*n+1)=DF_T_ij(mod((n)*n^2+n+2:n^2+1:n^4-(n)*n^2-2*n-n^2,n^2)); 

%% neumann boundary condition
A(n^2+2:n^2+1:(n-1)*n^2)=-3; % main diagonal 
A(2+(n+1)*n^2:n^2+1:(n^2)*(2*n-1))=4;  % mid
A(2+(2*n+1)*n^2:n^2+1:(n^2)*(3*n-1))=-1;  % mid

A(n^4-(n-2)*n^2-(n-2):n^2+1:n^4-n^2-1)=-3; % last squre matrix
A(n^4-(n+n-2)*n^2-(n-2):n^2+1:n^4-n^2-1-(n)*n^2)=4;  % last squre matrix
A(n^4-(2*n+n-2)*n^2-(n-2):n^2+1:n^4-n^2-1-(n*2)*n^2)=-1;  % last squre

%% dirichlet boundary condition  
A(1:(n^2)*(n)+n:n^4)=1;
A((n^2)*(n-1)+n:(n^2)*(n)+n:n^4)=1;

%% left diag Tij-1
A(n+2+n^2:n^2+1:n^4-(2*n)*n^2-n^2)=DF_T_ijminus1(mod((n+2+n^2:n^2+1:n^4-(2*n)*n^2-n^2),n^2));  
A((n-1)*(n^2)+n+n:(n^2)*(n)+n:n^4-(2*n-1)*n^2)=0;
A((n)*(n*n)+n+n+1:(n^2)*(n)+n:n^4-(2*n-1)*n^2)=0;

%% right diag Tij+1
A((n+n+1)*n^2+n+2:n^2+1:n^4-n-n^2)=DF_T_ijplus1(mod(((n+n+1)*n^2+n+2:n^2+1:n^4-n-n^2),n^2));  
A((2*n-1)*n^2+(n^3)+2*n:(n^2)*(n)+n:n^4-n-n^2)=0;
A((2*n)*n^2+(n^3)+2*n+1:(n^2)*(n)+n:n^4-n-n^2)=0;

%% left diag Ti-1j
A((n)*n^2+n+2:n^2+1:n^4-(n)*n^2-2*n-n^2)=DF_T_iminus1j(mod((n)*n^2+n+2:n^2+1:n^4-(n)*n^2-2*n-n^2,n^2)); 
A((n^3)+2*n+(n^2)*(n-2):(n^2)*(n)+n:n^4-(n)*n^2-2*n-n^2)=0;
A((n^3)+2*n+(n^2)*(n-1)+1:(n^2)*(n)+n:n^4-(n)*n^2-2*n-n^2)=0;

%% right diag Ti+1j
A((n+2)*n^2+n+2:n^2+1:n^4-(n)*n^2+n)=DF_T_iplus1j(mod(((n+2)*n^2+n+2:n^2+1:n^4-(n)*n^2+n),n^2)); % right mid
A((2*n)*n^2+n+n:(n^2)*(n)+n:n^4-(n)*n^2+n)=0;
A((2*n+1)*n^2+2*n+1:(n^2)*(n)+n:n^2*(n^2-n)+n)=0;

end

