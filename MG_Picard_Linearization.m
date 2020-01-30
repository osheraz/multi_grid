function [ one,two,three,four,five ] = MG_Picard_Linearization(dr,dQ,r,Q,dt,inx,n,one,two,three,four,five,K)
% function that calc the linerarized equation coefficient by Picard method
% recive grid information and K , return the seperated A matrix accordingly
% K can be changed in order to check diffrent problems.

K(:,:,:)=K(:,:,:).^2;
%K(:,:,:)=1;
             one(2:n(inx)-1,2:n(inx)-1,inx)= 3/(2*dt) + (K(3:n(inx),2:n(inx)-1,inx)+K(2:n(inx)-1,2:n(inx)-1,inx))./(2*dr(inx)^2)   +...
                                                 (K(2:n(inx)-1,2:n(inx)-1,inx)+K(1:n(inx)-2,2:n(inx)-1,inx))./(2*dr(inx)^2) +...
                                                 (K(2:n(inx)-1,3:n(inx),inx)+K(2:n(inx)-1,2:n(inx)-1,inx))./((2.*(dQ(inx)^2).*r(2:n(inx)-1,2:n(inx)-1,inx).^2)) +...
                                                 (K(2:n(inx)-1,2:n(inx)-1,inx)+K(2:n(inx)-1,1:n(inx)-2,inx))./((2.*(dQ(inx)^2).*r(2:n(inx)-1,2:n(inx)-1,inx).^2));
                                             
             two(2:n(inx)-1,2:n(inx)-1,inx)=    (K(2:n(inx)-1,2:n(inx)-1,inx))./(r(2:n(inx)-1,2:n(inx)-1,inx).*dr(inx)) +   (K(3:n(inx),2:n(inx)-1,inx)+K(2:n(inx)-1,2:n(inx)-1,inx))./(2*dr(inx)^2)   ;
           three(2:n(inx)-1,2:n(inx)-1,inx)=  - (K(2:n(inx)-1,2:n(inx)-1,inx))./(r(2:n(inx)-1,2:n(inx)-1,inx).*dr(inx)) +   (K(2:n(inx)-1,2:n(inx)-1,inx)+K(1:n(inx)-2,2:n(inx)-1,inx))./(2*dr(inx)^2)  ;
            four(2:n(inx)-1,2:n(inx)-1,inx)=    (K(2:n(inx)-1,2:n(inx)-1,inx))./((r(2:n(inx)-1,2:n(inx)-1,inx).^2).*(tan(Q(2:n(inx)-1,2:n(inx)-1,inx)).*(2*dQ(inx))))   +   (K(2:n(inx)-1,3:n(inx),inx)+K(2:n(inx)-1,2:n(inx)-1,inx))./((2.*(dQ(inx)^2).*r(2:n(inx)-1,2:n(inx)-1,inx).^2));
            five(2:n(inx)-1,2:n(inx)-1,inx)=  - (K(2:n(inx)-1,2:n(inx)-1,inx))./((r(2:n(inx)-1,2:n(inx)-1,inx).^2).*(tan(Q(2:n(inx)-1,2:n(inx)-1,inx)).*(2*dQ(inx))))   +   (K(2:n(inx)-1,2:n(inx)-1,inx)+K(2:n(inx)-1,1:n(inx)-2,inx))./((2.*(dQ(inx)^2).*r(2:n(inx)-1,2:n(inx)-1,inx).^2));
end

