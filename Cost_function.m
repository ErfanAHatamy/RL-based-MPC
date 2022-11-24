function Z = Cost_function(X0 , Ref , U , U0 , Ts , Np , Nc , A , B , C , D , gama , H , h , P, p)
Tp = Np * Ts;
time_pred = 0:Ts:Tp ;

X = zeros(size(X0,1) , numel(time_pred)-1) ;
Y = zeros(size([X0;U0],1) , numel(time_pred)-1) ;
X(: , 1) = X0 ;

for  i = 2:Np % It must be Np-1, however we begins the loop from 1 not 0, So it is true!
  X(:,i) = A * X(:,i-1) + B * U(i-1,:);
  Y(:,i) = [X(:,i);U(i-1,:)];
end

Error = reshape(Ref(:,2:Np) - Y(:,2:Np), [(Np-1)*size(Y,1),1]);
Error = [0;0;0;Error];  % The first 3 elements are X0 and U0 which are not variable!
% Z = Error' * H * Error + h' * Error + X(:,end)' * P * X(:,end) + p' * X(:,end)
Z = Error' * H * Error;

end


