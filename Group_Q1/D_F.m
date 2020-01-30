function [F_T,DF_T_ij,DF_T_iplus1j,DF_T_iminus1j,DF_T_ijplus1,DF_T_ijminus1] = D_F(Tn,n,r,Q,dr,dQ);
%D_F is a function which reply the equation F and its derivation by every
%parameter, after inserting the current T temperatures in it.
          
%define the parameters T_ij and its associates
T_ij_syms=Tn;
T_iplus1j_syms=circshift(Tn,-1);
T_iminus1j_syms=circshift(Tn,1);
T_ijplus1_syms=circshift(Tn,-n);
T_ijminus1_syms=circshift(Tn,n);
r_syms=r;
Q_syms=Q;


F_T=(- 2.*tan(Q_syms).*T_ij_syms.^2.*dQ.^2.*r_syms.^2 - 2.*tan(Q_syms).*T_ij_syms.^2.*dr.^2 + T_ij_syms.*T_ijplus1_syms.*dQ.*dr.^2 + 2.*tan(Q_syms).*T_ij_syms.*T_iplus1j_syms.*dQ.^2.*dr.*r_syms - T_ij_syms.*T_ijminus1_syms.*dQ.*dr.^2 - 2.*tan(Q_syms).*T_ij_syms.*T_iminus1j_syms.*dQ.^2.*dr.*r_syms + tan(Q_syms).*T_ijplus1_syms.^2.*dr.^2 + tan(Q_syms).*T_iplus1j_syms.^2.*dQ.^2.*r_syms.^2 + tan(Q_syms).*T_ijminus1_syms.^2.*dr.^2 + tan(Q_syms).*T_iminus1j_syms.^2.*dQ.^2.*r_syms.^2)./(2.*dQ.^2.*dr.^2.*r_syms.^2.*tan(Q_syms));
 
 
DF_T_ij=-(4.*T_ij_syms.*dr.^2.*tan(Q_syms) - T_ijplus1_syms.*dQ.*dr.^2 + T_ijminus1_syms.*dQ.*dr.^2 + 4.*T_ij_syms.*dQ.^2.*r_syms.^2.*tan(Q_syms) - 2.*T_iplus1j_syms.*dQ.^2.*dr.*r_syms.*tan(Q_syms) + 2.*T_iminus1j_syms.*dQ.^2.*dr.*r_syms.*tan(Q_syms))./(2.*dQ.^2.*dr.^2.*r_syms.^2.*tan(Q_syms));
 
 
DF_T_iplus1j=(T_ij_syms.*dr + T_iplus1j_syms.*r_syms)./(dr.^2.*r_syms);
 
 
DF_T_iminus1j=-(T_ij_syms.*dr - T_iminus1j_syms.*r_syms)./(dr.^2.*r_syms);
 
 
DF_T_ijplus1=(T_ij_syms.*dQ + 2.*T_ijplus1_syms.*tan(Q_syms))./(2.*dQ.^2.*r_syms.^2.*tan(Q_syms));
 
 
DF_T_ijminus1=-(T_ij_syms.*dQ - 2.*T_ijminus1_syms.*tan(Q_syms))./(2.*dQ.^2.*r_syms.^2.*tan(Q_syms));
 



end



