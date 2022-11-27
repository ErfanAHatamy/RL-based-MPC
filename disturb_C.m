function [C_ineq , C_eq] = disturb_C(teta, k, w)
global Nc Np
sx = size(teta.M,1);     % Size of state vector

C_in = zeros(Nc,sx);
for i = 1:k
    C_in(i,:) = (teta.M*w(:,i) - teta.m)';
end
C_ineq = C_in';

C_eq = [] ;

end