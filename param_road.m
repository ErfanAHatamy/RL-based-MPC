%% parameter list for MPC controller
% Author   : Shilp Dixit
% Date     : 01 November 2017
% Location : 17AA03 UoS

%% road properties
lw      = 3.5;      % lane width [m]
nlane   = 02;       % number of lane(s) [-]

rfront  = 100;      % sensor range front [m]
rrear   = -60;      % sensor range rear [m]

rho     = 1.225;    % air density [kg/m^3]

safePF  = 45;       % safe potential field limit [-]