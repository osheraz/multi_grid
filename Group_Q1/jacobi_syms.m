function [F_syms,DF_T_ij_syms,DF_T_iplus1j_syms,DF_T_iminus1j_syms,DF_T_ijplus1_syms,DF_T_ijminus1_syms,r_syms,Q_syms] = jacobi_syms(dr,dQ)
%   It takes the parameters Tn r teta etc. 
%   and gives back the values of jacobian derivatives depends on k value 
syms F(T_ij_syms,T_iplus1j_syms,T_iminus1j_syms,T_ijplus1_syms,T_ijminus1_syms,r_syms,Q_syms,dr,dQ)

% enter here F():
  F(T_ij_syms,T_iplus1j_syms,T_iminus1j_syms,T_ijplus1_syms,T_ijminus1_syms,r_syms,Q_syms,dr,dQ)=...
     (T_ij_syms./(r_syms*dr)+( (T_iplus1j_syms+ T_ij_syms))./(2* dr.^2)).* T_iplus1j_syms  + ...
     (-T_ij_syms./(r_syms*dr)+(( T_iminus1j_syms+ T_ij_syms))./(2* dr.^2)).* T_iminus1j_syms  +...
     (T_ij_syms./(r_syms.^2.*tan(Q_syms)*2.*dQ)+ ((T_ijplus1_syms+T_ij_syms))./(2*r_syms.^2.*dQ^2)).* T_ijplus1_syms +...
     (-T_ij_syms./(r_syms.^2.*tan(Q_syms)*2.*dQ)+ ((T_ijminus1_syms+T_ij_syms))./(2*r_syms.^2.*dQ^2)).* T_ijminus1_syms  +...
     (-((T_iplus1j_syms+ T_ij_syms))./(2* dr.^2)-( (T_iminus1j_syms+ T_ij_syms))./(2* dr.^2)- ((T_ijplus1_syms+T_ij_syms))./(2*r_syms.^2.*dQ^2)-((T_ijminus1_syms+T_ij_syms))./(2*r_syms.^2.*dQ^2)).* T_ij_syms;
 
% derivate the F by jacobian for the 5 syms parameters 

DF_T_ij_syms=simplify(jacobian(F,[T_ij_syms]));
DF_T_iplus1j_syms=simplify(jacobian(F,[T_iplus1j_syms]));
DF_T_iminus1j_syms=simplify(jacobian(F,[T_iminus1j_syms]));
DF_T_ijplus1_syms=simplify(jacobian(F,[T_ijplus1_syms]));
DF_T_ijminus1_syms=simplify(jacobian(F,[T_ijminus1_syms]));
F_syms=simplify(F);
r_syms=r_syms;
Q_syms=Q_syms;
end
