function [C_ineq , C_eq] = NLC(c_k, G, g, x, U, s_max2)
global C D Nc Np
sx = size(x,1);     % Size of state vector
C_in = zeros(Nc,size(C,2)+sx);
for i = 1:Np-1 % it should be Nc (instead of Np-1) when we use control horzion!
    C_in(i,1:sx) = (C*x(:,i) + D*U(i) + c_k(:,i))'; % it has traspose! 
    C_in(i,sx+1:end) = (x(:,i).^2 - s_max2)'; % it has traspose!
end
C_terminal = [(G*x(:,Nc)+ g)',(x(:,Np).^2 - s_max2)'];
C_ineq     = [C_in;C_terminal];

C_eq = [] ;

end