clear all
clc
nx=101;
n=nx-1;
dt=0.001;

rmin=0.5;
rmax=1;
Qmin=pi/6;
Qmax=5*pi/6;
r=linspace(rmin,rmax,n);
Q=linspace(Qmin,Qmax,n);
dr=r(2)-r(1); %const distance between points for r
dQ=Q(2)-Q(1); %const distance between point for Q
[Q,r]=meshgrid(Q(1:n),r(1:n));

%% init parameters
 load T0.mat
% Tn=linspace(0,1,n);% inital grid temp
% Tn=repmat(Tn',n,1);
% Tn=reshape(Tn,n,n);

%% init and sizing
T_old=Tn;
T_old_old=Tn;
one=zeros(n,n);
two=zeros(n,n);
three=zeros(n,n);
four=zeros(n,n);
five=zeros(n,n);
RHS=zeros(n,n);

%% creating the hole , updating DR,DQ
[DQ,DR,in,on,c,d,f,p] = creating_circle(dr,dQ,n,Q,r);
[Qf,rf]=pol2cart(Q,r); %change to polar

%% Coefficient calculation of the heat equation (with multiply dt)
[ one,two,three,four,five ] = Calc_Coefficient( n,DR,DQ,dt,dQ,dr,Q,r,one,two,three,four,five );

%% Time Loop
for t=0:dt:10*dt
    T_old_old(1:n,1:n)=T_old(1:n,1:n); % save fine grid n-2
    T_old(1:n,1:n)=  Tn(1:n,1:n); % save fine grid n-1
    RHS(1:end,1:end)= (-4.*T_old(1:end,1:end) + T_old_old(1:end,1:end)); %from jacobi
    norm=1;%init the norm
%% Iteration Loop
    while (norm>10^-4)
        T_check=Tn; % save fine grid k-1
        Tn(p)=1; Tn(f)=1 ;Tn(d)=1 ; Tn(c)=1; % initial the hole temp
        
        % neuman condition
        Tn(2:n-1,1)=(-Tn(2:n-1,3)+4.*Tn(2:n-1,2))./3;
        Tn(2:n-1,n)=(-Tn(2:n-1,n-3)+4.*Tn(2:n-1,n-2))./3;
        
        % jacobi iteration
        Tn(2:n-1,2:n-1)= ( (two(2:n-1,2:n-1)).* Tn(3:n,2:n-1)   +...
            (three(2:n-1,2:n-1)).* Tn(1:n-2,2:n-1)              +...
            (four(2:n-1,2:n-1)).* Tn(2:n-1,3:n)                 +...
            (five(2:n-1,2:n-1)).* Tn(2:n-1,1:n-2)               -...
            RHS(2:n-1,2:n-1) )./ one(2:n-1,2:n-1);
        Tn(in)=NaN;
        
        % update the convergence
        e=reshape(Tn-T_check,1,n^2);
        Tc=reshape(Tn,1,n^2);
        ind= find (abs(e)== max(abs(e)));
        norm = abs(e(ind))/abs(Tc(ind));
    end
end
%% plot fot evaluation
        contourf(Qf,rf,Tn(:,:),50);
        colorbar
        caxis([0,1])
        colormap jet
        