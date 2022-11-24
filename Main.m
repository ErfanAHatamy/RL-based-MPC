clc;clear all;close all;
%% ------------------------------------------------------------------------
% -------------------------- Initialization -------------------------------
T0 = 0;
Tf = 10;
Ts = 0.01;
t  = T0:Ts:Tf;
nt = numel(t);

% Dynamic Model parameters
A = [1 ,0.1;
     0 ,1];
B = [0.05; 0.1];
C_output = [1, 0;  
            0, 1]; % There may be different between state space output equation and inequality constarints, so two different matrices are defined.
D_output = [0; 0];
C = [0, 0;
     0, 0];
D = [0; 0];
c_bar = [0; 0];

rho = 0;
nu  = 0;


X = zeros(2,nt);
S = X;
X(:,1) = [rho; nu];
S(:,1) = [rho; nu];
u = zeros(1,nt);
a = u;

% MPC parameters
Np       = 20;
Nc       = Np;
rho_max  = 3; rho_min = -1;%!!!!!!!!!!!!!!!!!! modify limits 
nu_max   = 3; nu_min  = -1;%!!!!!!!!!!!!!!!!!! modify limits 
a_max    = 10; a_min  = -10;
s_max    = [rho_max; nu_max];
s_min    = [rho_min; nu_min];

H_states = [1;2];10*ones(size(X(:,1),1),1);
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
nu_r  = 1*ones(1,nt);
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
option = optimset('MaxIter',Opt_Iter,'Display','off','Algorithm',algorithm_name); % sqp,active-set,trust-region-reflective,interior-point

for i = 2:nt - Np
    c_k = c_bar + d_k;
    g_bar = c_k; % Remember that it can have more rows !! and there is a question about it: is it equal to c_bar or c_k?
    g_k = g_bar + h_k;
    NonLCon = @(U)NLC(C, D, c_k, G, g_k, X(:,i-1:i+Np-1), U, Nc, Np, s_max2) ; % Modify U size !!!!
    [Uopts , Fval(i), EXITFLAG , OUTPUT] = fmincon(@(U)Cost_function(X0 , Ref(:,i:i+Np-1), U , U0 , Ts , Np , Nc , A , B , C_output , D_output,  gama , H , h , P, p ) , Uopts , A_lin , B_lin , Aeq , Beq , LB , UB , NonLCon, option);

    u(1,i-1) = Uopts(1,1);
    U0 = Uopts(1,1);% Warm start
    X(:,i) = A * X(:,i-1) + B * u(:,i-1);
    X0 = X(:,i);
    number_of_iteration(i) = OUTPUT.iterations;
    
end


%% ------------------------------------------------------------------------
% ------------------------------- Plots -----------------------------------
figure(1);
subplot(3,1,1) ;
plot(t(1 , 1:nt-Np) , Ref(1 , 1:nt-Np) , t(1 , 1:nt-Np)  , X(1 , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('rho') ;
title('Trajectory Tracking') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')

subplot(3,1,2) ;
plot(t(1 , 1:nt-Np) , Ref(2 , 1:nt-Np) , t(1 , 1:nt-Np)  , X(2 , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('nu') ;
title('Trajectory Tracking') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')

subplot(3,1,3) ;
plot(t(1 , 1:nt-Np) , Ref(3 , 1:nt-Np) , t(1 , 1:nt-Np)  , u(: , 1:nt-Np) , 'LineWidth' , 1.25) ; hold on
xlabel('t (sec)') ;
ylabel('a') ;
title('Input Value') ;set(gca,'FontSize',10);
grid minor
legend('Reference Output','System Output')


figure(2);
plot(t(1:nt-Np)  ,Fval(1:nt-Np) , 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('Amp - CF') ;
title('Cost Function Variation') ;
grid on
legend('Cost');

figure(3);
plot(t(1:nt-Np)  ,number_of_iteration(1:nt-Np) , 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('number') ;
title('number of iteration') ;
grid on
legend('iterations');