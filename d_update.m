function d = d_update(theta,option_disturb)
global Np C Nk_prime
d = zeros(size(C,1), Np+1);

for k = 1:Np+1+Nk_prime 
    LCon = @(w)disturb_C(theta, k, w) ;
    w_opt = zeros(2,k);   
    [w1, d(1,k), ~ , ~] = fmincon(@(w)disturb(k, 1, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , d(2,k), ~ , ~] = fmincon(@(w)disturb(k, 2, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , d(3,k), ~ , ~] = fmincon(@(w)disturb(k, 3, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , d(4,k), ~ , ~] = fmincon(@(w)disturb(k, 4, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , d(5,k), ~ , ~] = fmincon(@(w)disturb(k, 5, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , d(6,k), ~ , ~] = fmincon(@(w)disturb(k, 6, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);

 
end


end
