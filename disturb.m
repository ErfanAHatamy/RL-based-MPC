function d = disturb(k , i, w)
global  A_cl C_cl
counter=zeros(2,1);
for j=1:k
    counter = counter + (A_cl^j)*w(:,j);
end
d = -1*C_cl(i,:)*counter;
end


