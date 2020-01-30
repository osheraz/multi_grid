clear all
close all;
%% problem parameters
nx=101;
n=nx-1;
rmin=0.5;
rmax=1;
Qmin=pi/6;
Qmax=5*pi/6;
r=linspace(rmin,rmax,n);
Q=linspace(Qmin,Qmax,n);
dr=r(2)-r(1);% calculation of generally dr
dQ=Q(2)-Q(1);% calculation of generally dQ
[Q,r]=meshgrid(Q(1:n),r(1:n));
[Qf,rf] = pol2cart(Q,r);

%% calc jacobain matrix with syms diff for build D_F function.
%[F_syms,DF_T_ij_syms,DF_T_iplus1j_syms,DF_T_iminus1j_syms,DF_T_ijplus1_syms,DF_T_ijminus1_syms,r_syms,Q_syms] = jacobi_syms_new(dr,dQ)

%% init parameters
A=sparse(n*n,n*n); % sparse Jacobian matrix init
r=reshape(r,n^2,1);
Q=reshape(Q,n^2,1);
Tn=linspace(0,1,n);% inital grid temp
Tn=repmat(Tn',n,1);
Tn=reshape(Tn,n^2,1);
F_T=0.*Tn;
DF_T_ij=0.*Tn;
DF_T_iplus1j=0.*Tn;
DF_T_iminus1j=0.*Tn;
DF_T_ijplus1=0.*Tn;
DF_T_ijminus1=0.*Tn;
dTn=0.*Tn;
norm=1;
i=0;
dTn_old=dTn;

%% plot inital temp distribution
Tn=reshape(Tn,n,n);
r=reshape(r,n,n);
subplot(2,2,1)
contourf(Qf,rf,Tn(:,:),50);
caxis([0,1])
colormap jet
title('Initial temperature distribution')
subplot(2,2,2)
plot(r(:,1),Tn(:,1));grid minor;xlabel('Radius[m]');ylabel('Temp[C]');
r=reshape(r,n^2,1);
Tn=reshape(Tn,n^2,1);

%% Iteration loop
while norm>1*10^-6
    
    T_old=Tn; %save Temp grid n-1
      
    %calc of jacobian matrix
    [F_T,DF_T_ij,DF_T_iplus1j,DF_T_iminus1j,DF_T_ijplus1,DF_T_ijminus1] = D_F(Tn,n,r,Q,dr,dQ);    
    %  create A matrix
    A=CreateMat(A,n ,DF_T_ij,DF_T_ijminus1,DF_T_ijplus1,DF_T_iminus1j,DF_T_iplus1j);
    %   create RHS
    RHS=CreateRHS(n,Tn,F_T);   
    % Normalize A,RHS
    [ A,RHS ] = Normalizing( A,RHS);
    % Solve Ax=b
    dTn=mldivide(A,-RHS);
    Tn=dTn+Tn;
    %   checking the convergence
    e=Tn'-T_old';
    ind= find (abs(e)== max(abs(e)));
    norm=abs(e(ind))/abs(Tn(ind)')
end
%% plot final temp distribution 
subplot(2,2,3)
Tn=reshape(Tn,n,n);
r=reshape(r,n,n);
Q=reshape(Q,n,n);
contourf(Qf,rf,Tn,40,'Linestyle','None');
title('final temperature distribution')
subplot(2,2,4)
plot(r(:,50),Tn(:,50));grid minor;xlabel('Radius[m]');ylabel('Temp[C]');

