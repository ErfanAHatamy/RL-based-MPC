function [Phi, Gamma, Psi, Omega] = func_predictionModel(A,B,Hp,Hc,Q,P,R)
    % func_predictionModel  Prediction Model for a Linear MPC controller
    %
    % [Phi, Gamma, Psi, Omega] = func_predictionModel(A,B,Hp,Hc,Q,P,R)
    %
    % Inputs
    %       A       : system matrix of a state-space system
    %       B       : input matrix of a state-space system
    %       Hp      : prediction horizon
    %       Hc      : control horizon
    %       Q       : stage state cost
    %       P       : terminal state cost
    %       R       : input penalty
    %
    % Outputs
    %       Phi     : stacked vector of system matrices
    %       Gamma   : stacked vector of input matrices
    %       Psi     : diagonal matrix of input penalty
    %       Omega   : diagonal matrix of state penalty
    %
    % Prediction Model:
    %       System Dynamics:
    %           x(k+1) = A*x(k) + B*u(k)
    % 
    %           [x(1);x(2);....;x(Hp)] = [A; A^2;...;A^(Hp)]*x(0) + [B, . . ,0;AB,B,...0;A^(Hp-1)*B,...]*Uk
    %           using x(0) = x(k)
    %       =>  Xk = Phi*x(k) + Gamma*Uk
    % 
    %       Cost Function:
    %           V(x,xt,u,zs) = (xi - xs)'*Q*(xi - xs) + (ui - us)'*R*(ui - us) +
    %                          (xN - xs)'*P*(xN-xs) + (xs - xt)'*T*(xs-xt)
    %           u and zs are decision variables
    %           x and xt are parameters
    % 
    %       Cost Function:
    %           (1/2)*u'*G*u + u'*F*x
    %           where
    %             u = [u;zs]
    %             x = [x,xt]
    %
    %	Author   : Shilp Dixit
    %   Date     : 27 October 2017
    %   Location : 17AA03 UoS

    %%
    nx = size(A,2);     % number of states
    nu = size(B,2);     % number of inputs

    % state prediction matrix
    Phi = zeros(Hp*nx, nx);
    Gamma1c = zeros(Hp*nx, nu);
    for ii = 1:Hp
        Phi((ii-1)*nx+1:ii*nx, 1:end) = A^(ii);

        Gamma1c((ii-1)*nx+1:ii*nx, 1:end) = A^(ii-1)*B;
    end

    % input prediction matrix
    Gamma = zeros(Hp*nx, Hc*nu);
    temp  = zeros(size(Gamma1c));
    for ii = 1:Hp
        temp = circshift(Gamma1c,(ii-1)*nx);
        if ii >= 2
            temp(1:(ii-1)*nx,:) = zeros((ii-1)*nx,nu);
        end

        if ii <= Hc
            Gamma(:,(ii-1)*nu+1:ii*nu) = temp;
        else
            Gamma(:,(Hc-1)*nu+1:Hc*nu) = Gamma(:,(Hc-1)*nu+1:Hc*nu) + temp;
        end
    end

    % input prediction matrix
    %{
    Gamma = zeros(Hp*nx, Hc*nu);
    for ii = 1:Hp
        if ii == 1
            Gamma((ii-1)*nx+1:ii*nx, 1:nu) = B;
        else
            flag = 1;
            for jj = 1:Hc
                if jj < Hc
                    Gamma((ii-1)*nx+1:ii*nx, (jj-1)*nu+1:jj*nu) = A^(ii-jj)*B;
                elseif jj >= Hc && flag == 1
                    temp = 0;
                    for kk = ii-Hc:-1:0
                        temp = A^(kk)*B + temp;
                    end
                    Gamma((ii-1)*nx+1:ii*nx, (jj-1)*nu+1:jj*nu) = temp;
                    flag = 0;
                end
            end
        end
    end
    %}

    Omega   = [];
    % stage and terminal cost (state)
    for ii = 1:Hp
        if ii < Hp
            Omega   = blkdiag(Omega, Q);
        else
            Omega   = blkdiag(Omega, P);
        end
    end

    % stage cost (input)
    Psi     = [];
    for ii = 1:Hc
        Psi = blkdiag(Psi, R);
    end

return