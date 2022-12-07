function [C_ineq , C_eq] = RL_C(M, m, w)
sw = size(w,2);
sm = size(M,1);
C_ineq = zeros(sw, sm);     % sw is variable but sm is 6

for i = 1:sw
    C_ineq(i,:) = M*w(:,i) - m;
end
C_eq =[];

end