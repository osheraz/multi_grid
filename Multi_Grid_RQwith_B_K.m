%Solution of heat equation  in a ring segment by a multigrid method
% Osher Azulay  , Assignment 1, Project 2 , Group 11.

clc;
clear all;
close all;

%% Change GridNumb,n in order to check others Grid 
GridNumb=5;   % Number of grids
n(1)=13;      % Coarse grid dimension
%% boundry elements
rmax=1;
rmin=0.5;
Qmax=5*pi/6;
Qmin=pi/6;
%% Create Grids
[ r,Q ,n,dr,dQ,Qf,rf] = MG_Create_Grids( GridNumb,rmin,rmax,Qmin,Qmax,n );
[Qf,rf] = pol2cart(Qf,rf);

%% init everything and sizing
Res(1:n(1),1:n(1),GridNumb)=0;
RHS(1:n(1),1:n(1),GridNumb)=0;
T_grid(1:n(1),1:n(1),GridNumb)=0;
T_grid_old(1:n(1),1:n(1),GridNumb)=0;
T_grid_old_old(1:n(1),1:n(1))=0;

%% init multipliers for every grid and Tc=0 for prolongation func
one(1:n(1),1:n(1),GridNumb)=0;
two(1:n(1),1:n(1),GridNumb)=0;
three(1:n(1),1:n(1),GridNumb)=0;
four(1:n(1),1:n(1),GridNumb)=0;
five(1:n(1),1:n(1),GridNumb)=0;
Tc(1:n(1),1:n(1),GridNumb)=0;

%% Change here initial state of fine grid 
%load T0_MG.mat  < with hole
%T_grid=T0_MG;
T_grid(1:n(1),1:n(1))=reshape(repmat(linspace(0,rmax,n(1)),1,n(1)),n(1),n(1)); % Linear
%T_grid(1,1:n(1),1)=0.1;
%T_grid(n(1),1:n(1),1)=0.1;

%% Initialization of the convergence normerion,iterations counter and old Grids
T_grid_old_old(1:n(1),1:n(1))=T_grid(1:n(1),1:n(1)); % save fine grid n-2
T_grid_old(1:n(1),1:n(1),1)=T_grid(1:n(1),1:n(1),1); % save fine grid n-1
dt=0.001;
k=3;
norm=1;
tic
counter=0;
%% Time Loop
for t=0:dt:10*dt  
    
    T_grid_old_old(1:n(1),1:n(1))=T_grid_old(1:n(1),1:n(1),1); % save fine grid n-2
    T_grid_old(1:n(1),1:n(1),1)=T_grid(1:n(1),1:n(1),1); % save fine grid n-1
    RHS(2:n(1)-1,2:n(1)-1,1)= ( -4.*T_grid_old(2:n(1)-1,2:n(1)-1,1) + T_grid_old_old(2:n(1)-1,2:n(1)-1))./(2*dt);
    norm=1;
%% V Cycles Loop   
    while norm>1e-4 
        
        T_grid(1:n(1),1:n(1),2:GridNumb)=0; % init the error in the coarse grids to 0.
        T_grid_check(1:n(1),1:n(1))=T_grid(1:n(1),1:n(1),1); % save fine grid k-1
        for inx=1:GridNumb
            % Calc Coefficient by Picard method
            [ one,two,three,four,five ] = MG_Picard_Linearization( dr,dQ,r,Q,dt,inx,n,one,two,three,four,five,T_grid_old);
            % Solve Linearized Picard equation by Jacobi method
            T_grid  = MG_jacobi(k,one,two,three,four,five,T_grid,RHS,n,inx,Qf,rf,dt);
            % Neuman boundry condtion
            T_grid(2:n(inx)-1,1,inx) = (4/3).*T_grid(2:n(inx)-1,2,inx)-(1/3).*T_grid(2:n(inx)-1,3,inx);
            T_grid(2:n(inx)-1,n(inx),inx) = (4/3).*T_grid(2:n(inx)-1,n(inx)-1,inx)-(1/3).*T_grid(2:n(inx)-1,n(inx)-2,inx);
            % calc of residual
            [ Res ] = MG_Calc_Res(Res,one,two,three,four,five,T_grid,RHS,n,inx);
            % calc R(r(inx-1)),R(T_grid_old(inx-1)) for all the coarse grids
            if (inx<GridNumb)
                RHS(1:n(inx+1),1:n(inx+1),inx+1)=Restriction( Res(1:n(inx),1:n(inx),inx),n(inx+1));
                T_grid_old(1:n(inx+1),1:n(inx+1),inx+1)=Restriction( T_grid_old(1:n(inx),1:n(inx),inx),n(inx+1));
            end
        end
        
        for inx=GridNumb-1:-1:1
            % calc T_grid with Prolongation of the error from T_grid_inx+1
            T_grid(1:n(inx),1:n(inx),inx)= T_grid(1:n(inx),1:n(inx),inx)+Prolongation(T_grid(1:n(inx+1),1:n(inx+1),inx+1),Tc(1:n(inx),1:n(inx),inx),n(inx),n(inx+1));
            T_grid  = MG_jacobi( k,one,two,three,four,five,T_grid,RHS,n,inx,Qf,rf,dt);
            % Neuman boundry condtion
            T_grid(2:n(inx)-1,1,inx) = (4/3).*T_grid(2:n(inx)-1,2,inx)-(1/3).*T_grid(2:n(inx)-1,3,inx);
            T_grid(2:n(inx)-1,n(inx),inx) = (4/3).*T_grid(2:n(inx)-1,n(inx)-1,inx)-(1/3).*T_grid(2:n(inx)-1,n(inx)-2,inx);
        end
 %% Checking the convergence   
 
        e=reshape((T_grid(1:n(1),1:n(1),1)-T_grid_check(1:n(1),1:n(1))),1,(n(1))^2);
        Tn=reshape(T_grid(1:n(1),1:n(1),1),1,(n(1))^2);
        ind= find (abs(e)== max(abs(e)));
        norm=abs(e(ind))/abs(Tn(ind))
        counter=counter+1;

    end
end
toc
disp(['V cycle counter: ' , int2str(counter)])

%% Plot for eval
contourf(Qf,rf,T_grid(:,:,1),50,'LineStyle','none');
caxis([0,1])
colorbar
colormap jet
pause(0.0001)
figure(2)
plot(r(:,20),T_grid(:,20))
xlabel('Radius[m]');ylabel('Temp[C]');
grid minor
