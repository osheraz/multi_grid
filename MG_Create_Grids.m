function [ r,Q,n,dr,dQ,Qf,rf] = MG_Create_Grids(GridNumb,rmin,rmax,Qmin,Qmax,n )
% function that creates the 3D Meshgrids from the user's input
% recive number of grid, boundry elements, n
% return r,Q,dr,dQ matrix ,sizes accordingly to users's input

%Coarse grid dimension
for i=2:GridNumb
    n(i)=2*n(i-1)-1;
end;

n=fliplr(n);
r(1:n(1),GridNumb)=0;
Q(1:n(1),GridNumb)=0;

%creating the coarse grid r and Q coordinates
r(1:n(GridNumb),GridNumb)=linspace(rmin,rmax,n(GridNumb)); 
Q(1:n(GridNumb),GridNumb)=linspace(Qmin,Qmax,n(GridNumb)); %creating the coarse grid Y coordinates

dr(GridNumb)=r(2,GridNumb)-r(1,GridNumb);
dQ(GridNumb)=Q(2,GridNumb)-Q(1,GridNumb);

for i=GridNumb:-1:1
    r(1:n(i),i)=linspace(rmin,rmax,n(i));
    Q(1:n(i),i)=linspace(Qmin,Qmax,n(i));
    dr(i)=r(2,i)-r(1,i);
    dQ(i)=Q(2,i)-Q(1,i);
    
end

[Qf,rf]=meshgrid(Q(:,1),r(:,1)); 

R(1:n(1),1:n(1),GridNumb)=0;
q(1:n(1),1:n(1),GridNumb)=0;
% init r and Q again with same val in r for colum, Q for row (3D Meshgrid)
for i=1:GridNumb
    R(1:n(i),1:n(i),i)=repmat(r(1:n(i),i),1,n(i));
    q(1:n(i),1:n(i),i)=repmat(Q(1:n(i),i),1,n(i))';
end

r=R;
Q=q;

end

