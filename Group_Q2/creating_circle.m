function [DQ,DR,in,on,c,d,f,p,Qv,rv] = creating_circle(dr,dQ,n,Q,r)
% function that create the distance Matrix  between each point in the grid from it
% neighbours DR,DQ with the hole conditions

% hole adjustment
r0 = 0.75;
Q0 = 4*pi/6;
R = 0.05;
DR=dr*ones(n,n); %matrix for distance between points (for r)
DQ=dQ*ones(n,n);

% Creating the circle hole geometry
L = linspace(0,2.*pi,110.*pi);
rv = (r0+R.*cos(L))';
Qv = (Q0+R.*sin(L))';

% Creating boolean grid of in/out the hole
[in,on] = inpolygon(Q,r,Qv,rv);
c=on; d=on; f=on; p=on;%just to initial size

for i=2:n-1 %loop for all the points 
    for j=2:n-1
        if ( in(i,j)==1 && in(i-1,j)==0) %insert to c the points (first line) above the hole 
            DR(i-1,j)=abs(r(i-1,j)-r0)-sqrt(R^2-(Q(i-1,j)-Q0).^2);
            c(i-1,j)=1;
        end
        if ( in(i,j)==1 && in(i+1,j)==0) %insert to d the points (first line) below the hole
            d(i+1,j)=1;
            DR(i+1,j)=abs(r(i+1,j)-r0)-sqrt(R^2-(Q(i+1,j)-Q0).^2);
        end
        if ( in(i,j)==1 && in(i,j+1)==0) %insert to f the points (first line) right the hole
            f(i,j+1)=in(i,j);
            DQ(i,j+1)=Q(i,j+1)-Q0-sqrt(R^2-abs(r(i,j+1)-r0).^2);
        end
        if ( in(i,j)==1 && in(i,j-1)==0) %insert to p the points (first line) left the hole
            p(i,j-1)=in(i,j);
            DQ(i,j-1)=Q0-Q(i,j-1)-sqrt(R^2-abs(r(i,j-1)-r0).^2);
        end
    end
end

DQ(in)=0; %initial zero in the hole 
DR(in)=0;

% shifting the points to be the first circle outside the hole
c=circshift(c,[+1 0]);
d=circshift(d,[-1 0]);
p=circshift(p,[0 +1]);
f=circshift(f,[0 -1]);

end

