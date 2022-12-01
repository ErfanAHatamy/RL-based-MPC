function q = q_update(theta,option_disturb)
global C Nk_prime
q = zeros(size(C,1),1+Nk_prime);

for k = 1:1+Nk_prime 
    LCon = @(w)disturb_C(theta, k, w) ;
    w_opt = zeros(2,k);   
    [w1q, q(1,k), ~ , ~] = fmincon(@(w)disturb_q(k, 1, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , q(2,k), ~ , ~] = fmincon(@(w)disturb_q(k, 2, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , q(3,k), ~ , ~] = fmincon(@(w)disturb_q(k, 3, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , q(4,k), ~ , ~] = fmincon(@(w)disturb_q(k, 4, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , q(5,k), ~ , ~] = fmincon(@(w)disturb_q(k, 5, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);
    [~ , q(6,k), ~ , ~] = fmincon(@(w)disturb_q(k, 6, w), w_opt, ...
        [],[],[],[],[],[], LCon,option_disturb);

 
end


end
