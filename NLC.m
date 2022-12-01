function [C_ineq , C_eq] = NLC(c_bar, d, q, x, U, s_max2)
global Ts Np Nc Nk_prime A B C D gama H h P p C_cl A_cl
sx = size(x,1);     % Size of state vector
C_in = zeros(Nc,size(C,1)+sx);
C_terminal = zeros(Nk_prime+1,size(C,1)+sx);

for i = 2:Np % it should be Nc (instead of Np) when we use control horzion! (Np not Np+1) because we have different terminal constarint! and it must count from 0 but we count from 1!
    c_k = c_bar -  d(:,i-1); % d multiplied on -1 because of negative value of its cost function(disturb)    
    C_in(i,1:size(C,1)) = (C*x(:,i-1) + D*U(i-1,1) + c_k - U(i-1,2:end)')'; % it has 2 traspose! (1: overal, 2: for U(i-1,2:end))
    C_in(i,size(C,1)+1:end) = ((A * x(:,i-1) + B * U(i-1,1)).^2 - s_max2)'; % it has traspose!   
end

% Terminal Constraints
for k_prime = Np:Np + Nk_prime % Np+1 means N in paper, becasue we must count from '0' bu we have to count from '1'
    c_k   = c_bar - 1*d(:,k_prime); % d multiplied on -1 because of negative value of its cost function(disturb)      
    g_bar = c_k;
    g_k   = g_bar - q(:,k_prime-Np+1); % q multiplied on -1 because of negative value of its cost function(disturb)
    G     = C_cl*A_cl^(k_prime-Np);
    C_terminal(k_prime-Np+1,1:size(C,1)) = (G*(A * x(:,Np) + B * U(Np,1)) + g_k - U(Np,2:end)')';% it has 2 traspose! (1: overal, 2: for U(i-1,2:end))
    C_terminal(k_prime-Np+1,size(C,1)+1:end) = ((A * x(:,Np) + B * U(Np,1)).^2 - s_max2)';
end

C_ineq     = [C_in;C_terminal];

C_eq = [] ;

end