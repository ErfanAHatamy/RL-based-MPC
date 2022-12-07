clc
i=5
NonLCon1 = @(U)NLC1(d, q, X, i , U, U_sigma_opt(1,:)) ;
[U_sigma_opt , Fval, ~ , out] = fmincon(@(U)Q1(S , Ref, 5, U, 0.1) , ...
                                                    U_sigma_opt(1:9,1:7) , A_lin , B_lin , Aeq , Beq , LB(1:9,:) , UB , NonLCon1, option)