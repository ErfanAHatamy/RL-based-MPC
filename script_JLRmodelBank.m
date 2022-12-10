%% script for creating a bank of kinematic bicycle models of a vehicle
% Author   : Shilp Dixit
% Date     : 14 February 2018
% Location : 17AA03 UoS

%% %%%%%%% ENVIRONMENT SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;
% clear;
close all;

%% %%%%%%% PLOT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot options
%plotOpts;

%% %%%%%%% CREATE PARAMETER ARRAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vp = (5/18).*[95:1:120];    % velocity array [m/s]

%% %%%%%%% INITIALIZE VEHICLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% road parameters
param_road;

% vehicle parameters
sv.M    = 2412.503;     % mass of vehicle [kg]
sv.L    = 4.78;         % length of vehicle [m]
sv.w    = 1.7;          % width of subject vehicle [m] (from IPG carmaker)
sv.lf   = 1.446;        % distance CoG to front axle [m]
sv.lr   = 1.477;        % distance CoG to rear axle [m]
sv.wb   = sv.lf+sv.lr;  % wheelbase [m] (from IPG carmaker)

sv.Cd   = 0.3581;       % coefficient of drag
sv.A    = 3.04;         % frontal area [m^2]

vnom    = mean(vp);     % nominal velocity [m/s]
vmax    = max(vp);      % maximum velocity [m/s]

Caf_nom = 1.088619810748564e+05;
Car_nom = 1.088619810748564e+05;

% discrete sampling time
ts      = 0.1;          % sampling time [s]

%% %%%%%%% LATERAL DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create model bank
lat.modelBank.A = cell(length(vp),1); % stack of system matrices
lat.modelBank.B = cell(length(vp),1); % stack of input matrices
lat.modelBank.C = cell(length(vp),1); % stack of ouput matrices
lat.modelBank.D = cell(length(vp),1); % stack of feedthrough matrices

C     = eye(2);
D     = zeros(2,1);
for ii = 1:length(vp)
    A = [0, vp(1,ii);...
         0, 0];

    B = [vp(1,ii)*sv.lr/sv.wb;...
         vp(1,ii)/sv.wb];

    csys = ss(A, B, C, D);          % CT state-space
    dsys = c2d(csys, ts, 'tustin'); % discrete time state-space

    % system matrices
    lat.modelBank.A{ii} = dsys.A;
    lat.modelBank.B{ii} = dsys.B;
    lat.modelBank.C{ii} = dsys.C;
    lat.modelBank.D{ii} = dsys.D;
end

%% calculate nominal model
[lat.Abar, lat.Bbar] = func_nominalModel(lat.modelBank);
lat.nx  = size(lat.Abar,1);
lat.nu  = size(lat.Bbar,2);

%% designing stable controller for nominal system
lat.Q = blkdiag(1e1, 1e2);    % state penalty
lat.R = blkdiag(1e1);         % input penalty

% compute control law by lqr
[lat.P,~,lat.K] = dare( lat.Abar,lat.Bbar,lat.Q,lat.R );    % terminal set penalty (P)
lat.p   = [0.60,0.1600];
    lat.p   = [0.72,0.1600];
lat.K   = place( lat.Abar,lat.Bbar,lat.p );
lat.Ak  = (lat.Abar-lat.Bbar*lat.K);                % system matrix for error dynamics

%% limits
% define state limits
lat.xmin = -[0.05*lw; 0.035];
lat.xmax =  [2*lw; 0.035];
    %
    lat.xlim    = [-lat.xmin; lat.xmax];
    lat.ax      = [-eye(lat.nx);eye(lat.nx)];
% polyhedron of state constraints
lat.X   = Polyhedron( 'A',lat.ax,'b',lat.xlim );

% define input limits
lat.umin = -[0.02];
lat.umax =  [0.02];
    %
    lat.ulim    = [-lat.umin; lat.umax];
    lat.au      = [-eye(lat.nu);eye(lat.nu)];
% polyhedron of input constraints
lat.U   = Polyhedron( 'A',lat.au,'b',lat.ulim );

% combine state and input constraints
lat.Z   = lat.X*lat.U;

%% obtain vertices of disturbance polyhedron
%%{
for ii = 1:length(vp)
    lat.W(ii) = (lat.modelBank.A{ii} - lat.Abar)*lat.X + (lat.modelBank.B{ii} - lat.Bbar)*lat.U;
end

% combine all disturbance sets
Vtotal = [];
for ii = 1:length(vp)
    Vtotal = vertcat(Vtotal,lat.W(ii).V);
end
lat.W = Polyhedron(Vtotal);

% plot disturbance set
%{
figure;
    h_Wplot = plot(lat.W);
        alpha(0.5);
        hold on;
%}

% create robust positively invariant set for error
loopCount   = 1e2;
epsilon     = 1e-3;
[h_Fsplot, lat.Fs, h_Fasplot, lat.Fas] = func_maxRPISet( lat.W, lat.Ak, loopCount, epsilon );
%     legend([h_Wplot(end) h_Fsplot h_Fasplot],'W','F_{s}','F(\alpha,s)');

    clear xlim
%{
figure
    plot(lat.W,'color','k', 'LineWidth',1); alpha(0.9);
        legend('$W$','Interpreter','Latex');
        xlabel('Y [m]');
        ylabel('\psi [rad.]');
        xlim([-0.3 0.3]);
        ylim([-0.09 0.09]);
        hold on;
    plot(lat.Fas,'color','y','LineWidth',1); alpha(0.4);
        xlabel('Y [m]');
        ylabel('\psi [rad.]');
        legend('$W$','$Z$','Interpreter','Latex');
%}

%%
% tightenend state constraints of nominal system
lat.Xbar = lat.X - lat.Fas;         	% tightened state set
lat.Ubar = lat.U - (-lat.K)*lat.Fas;   	% tightened input set
lat.Zbar = lat.Xbar*lat.Ubar;

%% %%%%%%% LONGITUDINAL DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create model bank
long.modelBank.A = cell(length(vp),1); % stack of system matrices
long.modelBank.B = cell(length(vp),1); % stack of input matrices
long.modelBank.C = cell(length(vp),1); % stack of ouput matrices
long.modelBank.D = cell(length(vp),1); % stack of feedthrough matrices

C     = eye(1);
D     = zeros(1,1);
A1 = 1e-3*(-1/2)*( rho*sv.A*sv.Cd/sv.M )*min(vp);
A2 = 1e-3*(-1/2)*( rho*sv.A*sv.Cd/sv.M )*max(vp);
for ii = 1:length(vp)
    theta1 = ( max(vp)-vp(1,ii) )/( max(vp)-min(vp) );
    theta2 = ( vp(1,ii)-min(vp) )/( max(vp)-min(vp) );

    A = A1*theta1 + A2*theta2;

    B = 1;

    csys = ss(A, B, C, D);          % CT state-space
    dsys = c2d(csys, ts, 'tustin'); % discrete time state-space

    % system matrices
    long.modelBank.A{ii} = dsys.A;
    long.modelBank.B{ii} = dsys.B;
    long.modelBank.C{ii} = dsys.C;
    long.modelBank.D{ii} = dsys.D;
end

%% calculate nominal model
[long.Abar, long.Bbar] = func_nominalModel(long.modelBank);
long.nx  = size(long.Abar,1);
long.nu  = size(long.Bbar,2);

%% designing stable controller for nominal system
long.Q = blkdiag(1e1);	% state penalty
long.R = blkdiag(1e1);  % input penalty

% compute control law by lqr
[long.P,~,long.K] = dare( long.Abar,long.Bbar,long.Q,long.R );  % terminal set penalty (P)
long.p  = 0.9049;
long.K  = place( long.Abar,long.Bbar,long.p );
long.Ak = (long.Abar-long.Bbar*long.K);               % system matrix for error dynamics

%% limits
% define state limits
long.xmin = -( -min(vp) );
long.xmax = max(vp);
    %
    long.xlim    = [-long.xmin; long.xmax];
    long.ax      = [-eye(long.nx);eye(long.nx)];
% polyhedron of state constraints
long.X   = Polyhedron( 'A',long.ax,'b',long.xlim );

% define input limits
long.umin = -[1.5];
long.umax =  [1.5];
    %
    long.ulim    = [-long.umin; long.umax];
    long.au      = [-eye(long.nu);eye(long.nu)];
% polyhedron of input constraints
long.U   = Polyhedron( 'A',long.au,'b',long.ulim );

% combine state and input constraints
long.Z   = long.X*long.U;

%% obtain vertices of disturbance polyhedron
% %{
for ii = 1:length(vp)
    long.W(ii) = (long.modelBank.A{ii} - long.Abar)*long.X + (long.modelBank.B{ii} - long.Bbar)*long.U;
end

% combine all disturbance sets
Vtotal = [];
for ii = 1:length(vp)
    Vtotal = vertcat(Vtotal,long.W(ii).V);
end
long.W = Polyhedron(Vtotal);

% plot disturbance set
%%{
figure;
    h_Wplot = plot(long.W);
        alpha(0.5);
        hold on;
%}

% create robust positively invariant set for error
loopCount   = 2e2;
epsilon     = 1e-5;
[h_Fsplot, long.Fs, h_Fasplot, long.Fas] = func_maxRPISet( long.W, long.Ak, loopCount, epsilon );
%     legend([h_Wplot(end) h_Fsplot h_Fasplot],'W','F_{s}','F(\alpha,s)');
% close all

%% %%%%%%% COMBINED SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine lateral and longitudinal systems (state-space model, disturbance, error, rpi sets)
nx = lat.nx + long.nx;      % net state dimension
nu = lat.nu + long.nu;      % net input dimension

% full model
Abar = blkdiag(lat.Abar,long.Abar);             % system matrix
Bbar = [zeros(lat.nx,long.nu), lat.Bbar;...     % input matrix
        long.Bbar, zeros(long.nx,long.nu)];

% control parameters
Q_e	= blkdiag( lat.Q,long.Q );
R_e = blkdiag( lat.R,long.R );
P_e = blkdiag( lat.P,long.P );
K_e = flipud( blkdiag( lat.K,long.K ) );	% nominal controller
Ak_e= Abar - Bbar*K_e;                      % closed-loop dynamics

% input and state-space
xmin = [lat.xmin;long.xmin];
xmax = [lat.xmax;long.xmax];
xlim = [-xmin;xmax];
ax	 = [-eye(nx);eye(nx)];
X    = Polyhedron( 'A',ax,'b',xlim );

umin = [long.umin;lat.umin];
umax = [long.umax;lat.umax];
ulim = [-umin;umax];
au	 = [-eye(nu);eye(nu)];
U    = Polyhedron( 'A',au,'b',ulim );

Z = X*U;

% net rpi error set
Fas = lat.Fas*long.Fas;

% tightened input and state-space
Xbar = X - Fas;         % tightened state space
Ubar = U - (-K_e)*Fas;  % tightened input space
Zbar = Xbar*Ubar;

% tracking set
xtmin = [-1;-1;-1].*Xbar.b(1:nx,:);
xtmax = [1;1;1].*Xbar.b(nx+1:end,:);
    xtmin = [-1;-1;-1].*X.b(1:nx,:);
    xtmax = [1;1;1].*X.b(nx+1:end,:);

xbarlim = Xbar.b;
ubarlim = Ubar.b;

return

%% %%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %{

figure
    subplot(2,1,1);
        plot(X); hold on; alpha(0.5);
        plot(Xbar,'color','c'); alpha(0.3);
        plot(Fas,'color','k');
            xlabel('Y [m]');
            ylabel('\psi [rad.]');
            zlabel('v_{x} [m/s]');
    subplot(2,1,2);
        plot(U); hold on; alpha(0.5);
        plot(Ubar,'color','c'); alpha(0.3);
            xlabel('a_{x} [m/s^2]');
            ylabel('\delta [rad.]');
%}