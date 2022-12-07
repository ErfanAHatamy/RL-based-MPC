function q = Q1(S , Ref , i , U, U0)
global  Ts Np A B C H P P_s
Tp = Np * Ts;
time_pred = 0:Ts:Tp ;
X0 = S(:,i-1);
X = zeros(size(X0,1) , numel(time_pred)) ;
Y = zeros(size([X0;0],1) , numel(time_pred)) ;
Ref = Ref(:,i:i+Np);
cx = size(C,1);
cy = size(Y,1);
X(:, 1) = X0 ;
X(:, 2) = A * X(:,1) + B * U0 ;
Y(:, 2) = [X(:, 2); U0];

for  i = 3:Np+1
  X(:,i) = A * X(:,i-1) + B * U(i-2,1);
  Y(:,i) = [X(:,i);U(i-2,1)];
end

Xn = X(:,Np+1);
Error = reshape(Ref(:,1:Np) - Y(:,1:Np), [cy*Np,1]);
UU    = reshape(U(:,2:end),[cx*(Np-1),1]);
q = Error' * H * Error + Xn' * P * Xn +  P_s(1,1:cx*(Np-1)) * UU  ;

end