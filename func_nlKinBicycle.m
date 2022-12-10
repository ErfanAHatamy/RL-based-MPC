function xdot = func_nlKinBicycle(t, x, time, ax, delta)
    % func_nlKinBicycle  simulate lateral and yaw dynamics of a vehicle
    % with a non-linear kinematic bicycle model
    %
    % xdot = func_nlKinBicycle(t, x, time, ax, delta)
    %
    % Inputs
    %       t       : system matrix of a state-space system
    %       x       : input matrix of a state-space system
    %       time    : prediction horizon
    %       ax      : control horizon
    %       delta   : stage state cost
    %
    % Outputs
    %       xdot    : stacked vector of system matrices
    %
    %	Author   : Shilp Dixit
    %   Date     : 14 February 2018
    %   Location : 17AA03 UoS

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v.lf = 1.446;      % distance of front axle from CoG [m]
    v.lr = 1.477;       % distance of rear axle from CoG [m]
    v.wb = v.lf+v.lr;   % wheelbase [m]
    v.Cd = 0.3581;      % coefficient of drag [-]
    v.A  = 3.04;        % frontal area [m^2]
    v.M  = 2412.503;    % vehicle mass [kg]

    rho  = 1.2250;      % density of air [kg/m^3]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE SPACE EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % interpolate input signals
    a_xi      = interp1(time, ax, t);    % Interpolate acc at times t
    delta_i   = interp1(time, delta, t); % Interpolate delta at times t

    % vehicle side-slip calculation
    beta = atan( (v.lr/(v.wb))*tan(delta_i) );     % vehicle side-slip angle [rad.]

    % vehicle lateral acceleration
%     ay  = ( x(4)^2/v.lr )*sin(beta);               % vehicle lateral acceleration [m/s^2]

    % (Schildbach - SCMPC for autonomous lane-change)
    xdot = [x(4)*cos( x(3)+beta );...
            x(4)*sin( x(3)+beta );...
          ((x(4)*cos(beta))/(v.wb))*tan(delta_i);...
%             x(4)/(v.lr)*sin(beta);...
            a_xi - 0*(1/2)*rho*v.Cd*v.A*( x(4)^2 )/v.M];