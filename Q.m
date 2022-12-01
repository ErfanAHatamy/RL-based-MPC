function q = Q(X0 , Ref , U , U0)
global  Ts  Np Nc A B C D gama H h P p P_s
Tp = Np * Ts;
time_pred = 0:Ts:Tp ;

X = zeros(size(X0,1) , numel(time_pred)) ;
Y = zeros(size([X0;U0],1) , numel(time_pred)) ;
X(: , 1) = X0 ;

for  i = 2:Np+1 % It must be Np-1, however we begins the loop from 1 not 0, So it is true!
  X(:,i) = A * X(:,i-1) + B * U(i-1,1);
  Y(:,i) = [X(:,i);U(i-1,1)];
end

Xn = X(:,Np+1);
Error = reshape(Ref(:,2:Np) - Y(:,2:Np), [(Np-1)*size(Y,1),1]);
UU    = reshape(U(:,2:end),[size(C,1)*Np,1]);
Error = [0;0;0;Error];  % The first 3 elements are X0 and U0 which are not variable!
% Z = Error' * H * Error + h' * Error + X(:,end)' * P * X(:,end) + p' * X(:,end)
q = Error' * H * Error + Xn' * P * Xn +  P_s * UU  ;

end


