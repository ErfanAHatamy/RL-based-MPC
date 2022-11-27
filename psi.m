function Z = psi(X, X0 , Ref , U , U0, teta)
global  Ts  Np Nc A B C D gama H h P p
l = loss_rl(X0, Ref((:,2:Np)), U0, Ref((3:Np)));
NonLCon = @(U)NLC(C, D, c_k, G, g_k, X, U, Nc, Np, s_max2) ;
[~ , v, ~ , ~] = fmincon(@(U)V(X0 , Ref, U , ...
                             Ts , Np , Nc , A , B , C_output , D_output,  gama , H , h , P, p ) , ...
                             [] , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);

[~ , q, ~ , ~] = fmincon(@(U)Q(X0 , Ref, U , ...
                               U0 , Ts , Np , Nc , A , B , C_output , D_output,  gama , H , h , P, p ) , ...
                               [] , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);

Z = (l + v - q)^2;
end


