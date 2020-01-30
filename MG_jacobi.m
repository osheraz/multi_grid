function [ T_grid ] = MG_jacobi( k,one,two,three,four,five,T_grid,RHS,n,inx,Qf,rf,dt)
%function that preform Jacobi iteration k times on the current linear
%problem Ax=b, recive the seperated A matrix (one:five) and RHS vector
%return the corrected T_grid(:,:,inx)

            for kk=1:k %Jacobi iteration for calculation of corrections on every grid
                T_grid(2:n(inx)-1,2:n(inx)-1,inx)=  ( (two(2:n(inx)-1,2:n(inx)-1,inx)).* T_grid(3:n(inx),2:n(inx)-1,inx)     +...
                                                    (three(2:n(inx)-1,2:n(inx)-1,inx)).* T_grid(1:n(inx)-2,2:n(inx)-1,inx)   +...
                                                     (four(2:n(inx)-1,2:n(inx)-1,inx)).* T_grid(2:n(inx)-1,3:n(inx),inx)     +...
                                                     (five(2:n(inx)-1,2:n(inx)-1,inx)).* T_grid(2:n(inx)-1,1:n(inx)-2,inx)   -...
                                                       RHS(2:n(inx)-1,2:n(inx)-1,inx) )./ one(2:n(inx)-1,2:n(inx)-1,inx);  
              %pause(dt)
              % contourf(Qf,rf,T_grid(:,:,1))
            end

end

