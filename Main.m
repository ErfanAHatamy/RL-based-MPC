clc;clear all;close all;
global  Ts Np Nc Nk_prime A B C D gama H h P p P_s A_cl C_cl
%% ------------------------------------------------------------------------
% -------------------------- Initialization -------------------------------
T0 = 0;
Tf = 10;
Ts = 0.1;
t  = T0:Ts:Tf;
nt = numel(t);

% Dynamic Model parameters
A = [1 ,0.1;
     0 ,1];
B = [0.05; 0.1];
C_output = [1, 0;  
            0, 1];
D_output = [0; 0];
C = [-1, 0;
     0, -1;
     1,  0;
     0,  1
     0,  0;
     0,  0];
D = [0; 0; 0; 0; 1; -1];
c_bar = [-1; -1; -1; -1; -10; -10];

rho = 0;
nu  = 0;
desired_poles = [-0.1,0.1];
K = place(A,B,desired_poles);
A_cl = (A - B*K);
C_cl = (C - D*K);
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
M_theta = [0,1;
          1,0;
          0,-1;
          -1,0];
m_theta = 1.5*[0.05;0.05;0.05;0.05];
theta.M = M_theta;
theta.m = m_theta;

X      = zeros(2,nt);
S      = X;
X(:,1) = [rho; nu];
S(:,1) = [rho; nu];
u      = zeros(1,nt);
sigma  = zeros(size(C,1),nt);
a      = u;

% MPC parameters
Np       = 30;
Nk_prime = 1; % it is for forming the terminal constraint which is used k_prime in eq(26),
% however by value=0 it means that we only have constrait for k=Np, no more!
Nc        = Np;
rho_max   = 1.5; rho_min = -1;
nu_max    = 1.5; nu_min  = -1;
a_max     = 10; a_min  = -10;
s_max     = [rho_max; nu_max];
s_min     = [rho_min; nu_min];
sigma_min = 0*ones(1,size(C,1)); % sigma must be greater than 0!
sigma_max = 10*ones(1,size(C,1));


H_states = [2;2];10*ones(size(X(:,1),1),1);
H_input  = 0*ones(size(u(:,1),1),1);
h_states = 0*ones(size(X(:,1),1),1);
h_input  = 0*ones(size(u(:,1),1),1);
P        = 0*eye(size(X(:,1),1));
p        = 0*ones(size(X(:,1),1),1);
P_sigma  = 1*ones(1,size(C,1));
gama     = 1;

A_lin = [];
B_lin = [];
Aeq = [];
Beq = [];
Opt_Iter = 30;
algorithm_name = 'active-set';


%% ------------------------------------------------------------------------
% ------------------------ Refrence Trajectory ----------------------------
rho_r = 1*((t>=25)&(t<=120)) + -1*((t<25)|(t>120));
% rho_r = 1*ones(1,nt);
nu_r  = 0*ones(1,nt);
a_r   = 0*ones(1,nt); % (?)

Ref = [rho_r;nu_r;a_r];

%% ------------------------------------------------------------------------
% -------------------------- Noise production -----------------------------
noise = zeros(size(X,1),nt-Np);
for i = 1:nt - Np
    noise(:,i) = cprnd(1,M,m)'; % It has transpose!
end

%% ------------------------------------------------------------------------
% --------------------------- MPC Controller ------------------------------
X0  = X(:,1);
a0  = a(:,1);
U0  = a0;
U_opts = 0 * ones(Nc , size(u,1)) ;
sigma_opt = 0 * ones(Np,size(C,1)); % remember that it includes sigma_n too.
U_sigma_opt = [U_opts,sigma_opt];

HH  = [H_states;H_input];
H   = diag(repmat(HH,Np,1));

hh   = [h_states;h_input];
h    = repmat(hh,Np,1);

P_s = reshape(repmat(P_sigma,Np,1),[1,size(C,1)*Np]);

LB = repmat([a_min,sigma_min] , Nc,1) ;  % lower bond of input signal
UB = [];repmat([a_max,sigma_max] , Nc,1) ;  % upper bond of input signal 
s_max2 = s_max.^2;

Fval = zeros(1,nt);
number_of_iteration = zeros(1,nt);
option = optimset('MaxIter',Opt_Iter,'Display','off','Algorithm',algorithm_name);
option_disturb = optimset('MaxIter',Opt_Iter,'Display','off','Algorithm',algorithm_name);

for i = 2:nt - Np
    if i == 2
        d = d_update(theta,option_disturb);
        q = q_update(theta,option_disturb);
    end
    
    NonLCon = @(U)NLC(c_bar, d, q, X(:,i-1:i+Np), U, s_max2) ;
    [U_sigma_opt , Fval(i), EXITFLAG , OUTPUT] = fmincon(@(U)Q(X0 , Ref(:,i:i+Np), U , U0) , ...
                                                    U_sigma_opt , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);
   
    u(1,i-1) = U_sigma_opt(1,1);
    sigma(:,i-1) = U_sigma_opt(1,2:end);
    error = (S(:,i-1) - X(:,i-1));
    X(:,i) = A * X(:,i-1) + B * u(:,i-1);
    a(:,i-1) = u(:,i-1)- K*error;
    U0 = a(:,i-1);
    S(:,i) = A * S(:,i-1) + B * a(:,i-1) + 0*noise(:,i);
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
%---------------------------------
figure(4);
plot(t(1:nt-Np)  ,sigma(:,1:nt-Np) , 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('sigma value') ;
title('sigma value for every step') ;
grid on
legend('sigma 1','sigma 2','sigma 3','sigma 4','sigma 5','sigma 6');