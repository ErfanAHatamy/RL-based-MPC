clc;clear all;close all;
global  Ts Np Nc A B C D gama H h P p A_cl C_cl
%% ------------------------------------------------------------------------
% -------------------------- Initialization -------------------------------
T0 = 0;
Tf = 3;
Ts = 0.01;
t  = T0:Ts:Tf;
nt = numel(t);

% Dynamic Model parameters
A = [1 ,0.1;
     0 ,1];
B = [0.05; 0.1];
C_output = [1, 0;  
            0, 1];
D_output = [0; 0];
C = [0, 0;
     0, 0];
D = [0; 0];
c_bar = [0; 0];

rho = 0;
nu  = 0;
desired_poles = [-0.1,0.1];
K = place(A,B,desired_poles);
A_cl = (A - B*K);
C_cl = (C_output - D_output*K);
% Set Real W: Mw<m
M = [0,1;
     1,1;
     1,0;
     1,-1;
     0,-1;
     -1,-1;
     -1,0;
     -1,1];
m = [0.01;0.015;0.01;0.015;0.01;0.015;0.01;0.015];
M_teta = [0,1;
          1,0;
          0,-1;
          -1,0];
m_teta = 10*[0.05;0.05;0.05;0.05];
teta.M = M_teta;
teta.m = m_teta;

X = zeros(2,nt);
S = X;
X(:,1) = [rho; nu];
S(:,1) = [rho; nu];
u = zeros(1,nt);
a = u;

% MPC parameters
Np       = 20;
Nc       = Np;
rho_max  = 3; rho_min = -1;
nu_max   = 3; nu_min  = -1;
a_max    = 10; a_min  = -10;
s_max    = [rho_max; nu_max];
s_min    = [rho_min; nu_min];

H_states = [20;2];10*ones(size(X(:,1),1),1);
H_input  = 0*ones(size(u(:,1),1),1);
h_states = 0*ones(size(X(:,1),1),1);
h_input  = 0*ones(size(u(:,1),1),1);
P        = 0*eye(size(X(:,1),1));
p        = 0*ones(size(X(:,1),1),1);
gama     = 1;

A_lin = [];
B_lin = [];
Aeq = [];
Beq = [];
Opt_Iter = 30;
algorithm_name = 'active-set';


%% ------------------------------------------------------------------------
% ------------------------ Refrence Trajectory ----------------------------
% rho_r = 1*((t>=25)&(t<=120)) + -1*((t<25)|(t>120));
rho_r = 1*ones(1,nt);
nu_r  = 0*ones(1,nt);
a_r   = 0*ones(1,nt); % (?)

Ref = [rho_r;nu_r;a_r];

%% ------------------------------------------------------------------------
% --------------------------- MPC Controller ------------------------------
X0  = X(:,1);
a0  = a(:,1);
U0  = a0;
Uopts = zeros(Nc , size(u,1)) ;

% HH  = [H_states, zeros(size(H_states,1),size(H_input,2));
%        zeros(size(H_input,1),size(H_states,2)),H_input];
HH  = [H_states;H_input];
H   = diag(repmat(HH,Np,1));

hh   = [h_states;h_input];
h    = repmat(hh,Np,1);
d = zeros(2, Np+1);
for k = 1:Np+1
    LCon = @(w)disturb_C(teta, k, w) ;
    best_guess = zeros(2,k);
    [~ , d(1,k), ~ , ~] = fmincon(@(w)disturb(k, 1, w), best_guess, ...
        [],[],[],[],[],[], LCon);
    [~ , d(2,k), ~ , ~] = fmincon(@(w)disturb(k, 2, w), best_guess, ...
        [],[],[],[],[],[], LCon);
end
d_k = [0; 0];
G   = [0, 0;
       0, 0];
g_k = [0; 0];
h_k = [0; 0];

LB = repmat(a_min , Nc,1) ;  % lower bond of input signal
UB = repmat(a_max , Nc,1) ;  % upper bond of input signal
s_max2 = s_max.^2;

Fval = zeros(1,nt);
number_of_iteration = zeros(1,nt);
option = optimset('MaxIter',Opt_Iter,'Display','off','Algorithm',algorithm_name);

for i = 2:nt - Np
    c_k = c_bar + d;
    g_bar = c_bar + d_k; 
    g_k = g_bar + h_k;
    NonLCon = @(U)NLC(c_k, G, g_k, X(:,i-1:i+Np-1), U, s_max2) ;
    [Uopts , Fval(i), EXITFLAG , OUTPUT] = fmincon(@(U)Q(X0 , Ref(:,i:i+Np-1), U , U0) , ...
                                                    Uopts , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);
    %[Uopts , Fval(i), EXITFLAG , OUTPUT] = fmincon(@(U)V(X0 , Ref(:,i:i+Np-1), U) , ...
    %                                                Uopts , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);
    
    u(1,i-1) = Uopts(1,1);
    error = (S(:,i-1) - X(:,i-1))
    X(:,i) = A * X(:,i-1) + B * u(:,i-1);
    a(:,i-1) = u(:,i-1)- K*error;
    U0 = a(:,i-1);
    S(:,i) = A * S(:,i-1) + B * a(:,i-1) + cprnd(1,M,m)';
    X0 = S(:,i);
    number_of_iteration(i) = OUTPUT.iterations;
    
end


%% ------------------------------------------------------------------------
% ------------------------------- Plots -----------------------------------
figure(1);
subplot(3,1,1) ;
plot(t(1 , 1:nt-Np) , Ref(1 , 1:nt-Np) , t(1 , 1:nt-Np)  , S(1 , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('rho') ;
title('Trajectory Tracking') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')

subplot(3,1,2) ;
plot(t(1 , 1:nt-Np) , Ref(2 , 1:nt-Np) , t(1 , 1:nt-Np)  , S(2 , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('nu') ;
title('Trajectory Tracking') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')

subplot(3,1,3) ;
plot(t(1 , 1:nt-Np) , Ref(3 , 1:nt-Np) , t(1 , 1:nt-Np)  , a(: , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('a') ;
title('Input Value') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')

%---------------------------------
figure(2);
plot(t(1:nt-Np)  ,Fval(1:nt-Np) , 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('Amp - CF') ;
title('Cost Function Variation') ;
grid on
legend('Cost');
%---------------------------------
figure(3);
plot(t(1:nt-Np)  ,number_of_iteration(1:nt-Np) , 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('number') ;
title('number of iteration') ;
grid on
legend('iterations');