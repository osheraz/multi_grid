function FgProlong = Prolongation(CgCorr,FgProlong,Nf,Nc)

%for i=3:2:Nf-2
%    for j=2:2:Nf-1
%        FgProlong(i,j)=(CgCorr((i-1)/2+1,j/2)+CgCorr((i-1)/2+1,j/2+1))/2;
%    end
%end

%%  Prolongation(T_grid(1:n(l+1),1:n(l+1),l+1),Tc(1:n(l),1:n(l),l),n(l),n(l+1));  

FgProlong(3:2:Nf-2,2:2:Nf-1)=(CgCorr(2:Nc-1,1:Nc-1)+CgCorr(2:Nc-1,2:Nc))/2;            
%for i= 2:2:Nf-1  
%    for j= 3:2:Nf-2
%        FgProlong(i,j)=(CgCorr( i/2,(j-1)/2+1)+CgCorr(i/2+1,(j-1)/2+1))/2;
%    end
%end
FgProlong(2:2:Nf-1,3:2:Nf-2)=(CgCorr(1:Nc-1,2:Nc-1)+CgCorr(2:Nc,2:Nc-1))/2;
%for i= 2:2:Nf-1  
%    for j= 2:2:Nf-1
%        FgProlong(i,j)=0.25*(CgCorr( i/2,j/2)+CgCorr( i/2+1,j/2)+CgCorr( i/2,j/2+1)+CgCorr( i/2+1,j/2+1));
%    end
%end
FgProlong(2:2:Nf-1,2:2:Nf-1)=0.25*(CgCorr(1:Nc-1,1:Nc-1)+CgCorr(2:Nc,1:Nc-1)+CgCorr(1:Nc-1,2:Nc)+CgCorr(2:Nc,2:Nc));
%for i=1:2:Nf
%    for j=1:2:Nf
%        FgProlong(i,j)=CgCorr((i+1)/2,(j+1)/2);
%    end
%end
FgProlong(3:2:Nf-2,3:2:Nf-2)=CgCorr(2:Nc-1,2:Nc-1);
end

