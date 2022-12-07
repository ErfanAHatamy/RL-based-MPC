function [C_ineq , C_eq] = NLC1(d, q, X, i, U, U0)
global Np Nk_prime A B C D C_cl A_cl c_bar s_max2
U1 = [[U0,zeros(1,6)];U];
sx = size(X,1);     % Size of state vector
sc = size(C,1);
C_in = zeros(Np,sc+sx);
C_terminal = zeros(Nk_prime+1,sc+sx);
x = X(:,i-1:i+Np-1);
for i = 1:Np
    c_k = c_bar -  d(:,i); % d multiplied on -1 because of negative value of its cost function(disturb)    
    C_in(i,1:sc) = (C*x(:,i) + D*U1(i,1) + c_k - U1(i,2:end)')'; % it has 2 traspose! (1: overal, 2: for U(i-1,2:end))
    C_in(i,sc+1:end) = ((A * x(:,i) + B * U1(i,1)).^2 - s_max2)'; % it has traspose!   
end

% Terminal Constraints
for k_prime = Np:Np + Nk_prime % Np+1 means N in paper, becasue we must count from '0' bu we have to count from '1'
    c_k   = c_bar - 1*d(:,k_prime); % d multiplied on -1 because of negative value of its cost function(disturb)      
    g_bar = c_k;
    g_k   = g_bar - q(:,k_prime-Np+1); % q multiplied on -1 because of negative value of its cost function(disturb)
    G     = C_cl*A_cl^(k_prime-Np);
    C_terminal(k_prime-Np+1,1:sc) = (G*(A * x(:,Np) + B * U1(Np,1)) + g_k - U1(Np,2:end)')';% it has 2 traspose! (1: overal, 2: for U(i-1,2:end))
    C_terminal(k_prime-Np+1,sc+1:end) = ((A * x(:,Np) + B * U1(Np,1)).^2 - s_max2)';
end

C_ineq     = [C_in;C_terminal];

C_eq = [] ;

end