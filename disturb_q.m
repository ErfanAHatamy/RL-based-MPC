function d = disturb_q(k , i, w)
global  A_cl C_cl
counter=zeros(2,1);
for j=1:k
    counter = counter + (A_cl^j)*w(:,j);
end
G = C_cl*A_cl^(k-1);
d = -1*G(i,:)*counter;
end


