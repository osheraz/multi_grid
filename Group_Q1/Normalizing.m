function [ A,RHS ] = Normalizing( A,RHS)
% Normalize A,RHS by diag(A) in order to lower cond(A) number.
RHS=RHS./diag(A);
A=bsxfun(@rdivide, A, diag(A));

end

