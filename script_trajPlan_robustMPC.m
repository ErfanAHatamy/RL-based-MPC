%% script robust trajectory planning leveraging V2X - MPC design
% Author   : Shilp Dixit
% Date     : 14 February 2018
% Location : 17AA03 UoS

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENVIRONMENT SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

% return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plotOpts;
% quadprog optimisation setting
options = optimset('Algorithm','interior-point-convex',...
                   'Display','off');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD IPG CARMAKER RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('sim_leadVehicle');
    a = simLoad;
% create time array
time = a.Time.data;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MPC CONTROLLER PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MPC parameters
Hp = 25;    % Prediction horizon [-]
Hc = 01;    % Control horizon [-]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE MODEL, CONSTRAINT SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_JLRmodelBank;

% load reachability set
load reachBndry;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE NOMINAL DISCRETE LINEAR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
x(k+1) = A x(k) + B u(k)
y(k)   = C x(k) + D u(k)
%}

% system matrices
A   = Abar;
B   = Bbar;

% state, input and output dimension
nx  = size( A,1 );
nu  = size( B,2 );
nw  = nx;
% standard basis vector in disturbance space
e   = eye(nw);

% steady-state conditions
%{
(xk,uk) belong to Z = {z = [xk;uk] : Az*z <= bz}
zs = {[xs;us] = Mt*theta : Mt*theta in Z}

zs = Mt*theta
yt = Nt*theta
%}
Mt = [1, 0, 0;...
      0, 1, 0;...
      0, 0, 1;...
      0, 1, 0;...
      0, 1, 0];
% Mt = eye(nx+nu);
Mx = Mt(1:nx, :);
Mu = Mt(nx+1:end, :);
np = size(Mt, 2);

%{
compute positive definite matrix P by soliving discrete lyapunov equation
(A+B*K)'*P*(A+B*K) - P = -(Q+K'*R*K)
%}
Q = blkdiag(2e-2, 1e-2, 1e1);       % state penalty
R = blkdiag(1.5e0, 2e2);            % input penalty

% compute control law by lqr
[P,~,K] = dare( A,B,Q,R );  % terminal set penalty (P)
T       = 1e2*P;            % tracking error

Ak      = (A-B*K);          % system matrix for error dynamics

% closed loop system
E = [A - eye(nx), B;...
     eye(nx)    , zeros(nx,nu)];

if rank(E) < nx + nu
    disp('free variables required');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE PREDICTION MODEL AND COST FUNCTION MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create prediction model matrices
[Phi, Gamma, Psi, Omega] = func_predictionModel( A,B,Hp,Hc,Q,P,R );
% [H,S] = create_CHS(A,B,eye(nx),Hp,Hc);

% return

Jx  = repmat( eye(nx), Hp, 1 );
Ju  = repmat( eye(nu), Hc, 1 );

% Initialise cost function matrices
    G11 = 2*( Q + Phi'*Omega*Phi );
    G12 = zeros( nx,Hc*nu );
    G13 = zeros( nx,np );
    G21 = 2*2*Gamma'*Omega*Phi;
    G22 = 2*(Psi + Gamma'*Omega*Gamma);
    G23 = 2*( -2*(Gamma'*Omega*Jx*Mx + Psi*Ju*Mu) );
    G31 = 2*( -2*Mx'*(Q + Jx'*Omega*Phi) );
    G32 = zeros( np,Hc*nu);
    G33 = 2*(Mt'*blkdiag(Hp*Q+P+T,Hp*R)*Mt);
G   = [G11, G12, G13;...
       G21, G22, G23;...
       G31, G32, G33];
G   = (G+G')/2;

    F11 = zeros( nx,nx );
    F12 = zeros( nx,nx );
    F21 = zeros( Hc*nu,nx );
    F22 = zeros( Hc*nu,nx );
    F31 = zeros( np,nx );
    F32 = -2*Mx'*T;
F   = [F11, F12;...
       F21, F22;...
       F31, F32];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE INVARIANT TERMINAL SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    defining Ktheta = [-K I]*Mt
    control law is defined by u = K*x + Ktheta*theta
%}
Ktheta = [K eye(nu)]*Mt;

% maximal invariant set for tracking
Ae = [ax     ,  zeros(size(ax,1),np);...
     -(A-B*K), -B*Ktheta;...
     +(A-B*K), +B*Ktheta;...
      K      , -Ktheta;...
     -K      ,  Ktheta;...
      zeros(2*(nx+nu),nx),[-Mx;Mx;-Mu;Mu]];

lambda = 0.99;
be = [Xbar.b;Xbar.b;Ubar.b;lambda*[Xbar.b;Ubar.b]];

Xfw = Polyhedron(Ae,be);                % maximal invariant set
Xf  = Xfw.projection(1:nx,'fourier');	% projection of maximal invariant set on system states
% return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESIGN STACKED CONSTRAINT MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constraint I (optimal initialisation of nominal system)
L1  = [-Fas.A, zeros(size(Fas.A,1),Hc*nu), zeros(size(Fas.A,1),np)];
c1  = Fas.b;
W1  = [-Fas.A, zeros(size(Fas.A,1),nx)];

% constraint II (state constraints of initial state of nominal system)
L2  = [Xbar.A, zeros(size(Xbar.A,1),Hc*nu), zeros(size(Xbar.A,1),np)];
c2  = Xbar.b;
W2  = zeros( size(Xbar.A,1),nx+nx );

% constraint III (state constraints of prediction model of nominal system - excluding terminal state)
L3  = [];
c3  = [];
W3  = [];
for ii = 1:size(Xbar.A,1)
    for jj = 1:nx:nx*(Hp-1)
        % extract relevant prediction matrices
        temp    = [ Phi(jj:jj+nx-1,:), Gamma(jj:jj+nx-1,:) ];

        trunL3  = zeros( 1,nx+(Hc*nu)+np );
        for kk = 1:nx
            trunL3 = trunL3 + [ Xbar.A(ii,kk).*temp(kk,:) zeros(1,np) ];
        end
        L3 = vertcat( L3, trunL3 );
        c3 = vertcat( c3,Xbar.b(ii,1) );
        W3 = vertcat( W3,zeros(size(Xbar.A(ii,:),1),nx+nx) );
    end
end

% constraint IV (input constraints for predicted inputs)
c4  = [];
tempL4 = [];
for ii = 1:Hc
    tempL4 = blkdiag( tempL4, Ubar.A );
    c4     = vertcat( c4, Ubar.b );
end
L4  = [ zeros(size(tempL4,1),nx), tempL4, zeros(size(tempL4,1),np) ];
W4  = zeros(size(tempL4,1),nx+nx);

% constraint V (terminal set constraint and steady-state constraint)
L5  = [];
c5  = [];
W5  = [];
for ii = 1:size(Xfw.A,1)
    L5 = vertcat( L5,[ Xfw.A(ii,1:nx)*[Phi(end-nx+1:end,:), Gamma(end-nx+1:end,:)], Xfw.A(ii,nx+1:end) ] );
    c5 = vertcat( c5,Xfw.b(ii,1) );
    W5 = vertcat( W5,zeros(size(Xfw.A(ii,:),1),nx+nx) );
end

% constraint VI (condition problem better by adding constraints for x_hat (reference))
L6  = [zeros(size(X.A,1),nx), zeros(size(X.A,1),Hc*nu), zeros(size(X.A,1),np)];
c6  = X.b;
W6  = [ zeros( size(X.A,1),nx ), -X.A ];

% combine all constraints
L   = vertcat(L1,L2,L3,L4,L5,L6);
c   = vertcat(c1,c2,c3,c4,c5,c6);
W   = vertcat(W1,W2,W3,W4,W5,W6);

% clear temp variables
clear temp trunL3 tempL3 tempL4

% return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAJECTORY PLANNING (AFP + MPC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial value
x0 = [a.Car_BdySensor_CG_Pos_0_y.data(1, 1) + lw/2;...
      a.Car_Yaw.data(1, 1);...
      0.98*max(vp)];

% creating simulation clock
t      = downsample(time, 2);
tSNAP  = min(t):ts:(max(t)+ts);
tictoc = zeros(1, length(t));

% initializing traffic position, velocity
lv.v = downsample(a.Traffic_LV_LongVel.data, 2); % velocity of lead vehicle [m/s]
lv.x = downsample(a.Traffic_LV_tx.data, 2);      % x-position of lead vehicle
    lv.x = lv.x - 0*115;
lv.y = downsample(a.Traffic_LV_ty.data, 2)+lw/2; % y-position of lead vehicle
lv.w = 1.7;     % width of lead vehicle [m]
lv.l = 4.1;     % length of lead vehicle [m]

h = Hp*ts;      % headway time [s]

% initialize SV state vector
x = [x0, zeros(nx, length(t)-1)];
s = zeros(1, length(t));	% longitudinal displacement [m]

% target states
xt = [zeros(1, length(t));...
      zeros(1, length(t));...
      vnom*ones(1, length(t))];

% initialize input, disturbance, and output vector
u  = zeros(nu, length(t));

clear xlim
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsnap   = [];
u_rob   = zeros(nu,length(t));
exitFlg = zeros(length(t),1);       % status of quadprog optimisation

for ii = 1:length(t)-1 % 125
    %  Run vehicle perception system

    % define origin of road co-ordinates
    ro.x = s(:, ii);
    ro.y = 0;

    % generate PF of immediate road ahead
    Road = func_getRoadPF(lw, nlane);

    % check if obstacle within range of SV
    if lv.x(1, ii) - s(:, ii) < rfront && lv.x(1, ii) - s(:, ii) > rrear
        lv.xpos = lv.x(1, ii) - s(:, ii);
        lv.ypos = lv.y(1, ii) + lw/2;

        % generate obstacle PF
        obs = func_getCarPF(x(3, ii),...
                         lv.v(1,ii),...
                         h,...
                         lv.l,...
                         lv.w,...
                         lv.xpos,...
                         lv.ypos,...
                         Road.Xg,...
                         Road.Yg);
    else
        lv.xpos = nan;
        lv.ypos = lv.y(1, ii) + lw/2;
        obs.pf = 0;
    end
    % exception handling
    if isempty(obs.pf)
        obs.pf = zeros(size(Road.pf));
    end
    % Combine obstacle pf with road pf to get overall PF
    Road.pf = Road.pf + obs.pf;
    % truncate very large values at road edges
    Road.pf(Road.pf > 150) = 150;
    Road.pf = round(Road.pf.*1e4)/1e4;

    % get SV current position (in road co-ordinates)
    sv.x = 0;
    sv.y = x(1, ii) + lw/2;
    % create subject vehicle outline
    sv.block = func_vehicleOutline( sv.x, sv.y, x(2,ii), sv.lf, sv.lr, sv.w );

    % LV basic shape for top-view (in road coordinates)
    lv.xpos = lv.x(1, ii) - s(:, ii);
    % create lead vehicle outline
    lv.block = func_vehicleOutline( lv.xpos, lv.ypos, 0, lv.l/2, lv.l/2, lv.w );

    % reachable set of cruising vehicle
    sv.reach.max.x = vnom*h;                            % maximum x-distance travelled [m]
    sv.reach.max.y = x(1, ii) + lw/2 + 3.992;           % left edge of road [m]
    sv.reach.min.x = vnom*h + (1/2)*umin(1, 1)*h^2;     % minimum x-distance travelled [m]
    sv.reach.min.y = max(lw/2, x(1, ii) + lw/2 - 3.992);% right edge of road [m]

    % extract indices of vertices of reachable set
    sv.reach.max.xi = find(Road.Xg(1, :) <= sv.reach.max.x, 1, 'last');
    sv.reach.max.yi = find(Road.Yg(:, 1) <= sv.reach.max.y, 1, 'last');
    sv.reach.min.xi = find(Road.Xg(1, :) <= sv.reach.min.x, 1, 'last');
    sv.reach.min.yi = find(Road.Yg(:, 1) <= sv.reach.min.y, 1, 'last');

    %%
    % generate PF of immediate road ahead
    Road.f = func_getRoadPF(lw, nlane);
    Road.f.pf(isinf(Road.f.pf)) = nan;
    % check if obstacle within range of SV
    if lv.x(1, ii) + lv.v(:, ii)*h*0.75 - s(:, ii) < rfront && lv.x(1, ii) + lv.v(:, ii)*h*0.75 - s(:, ii) > rrear
        lv.f.xpos = lv.x(1, ii) - s(:, ii) + lv.v(:, ii)*h*0.75;
        lv.f.ypos = lv.y(1, ii) + lw/2;

        % generate obstacle PF
        obs.f  = func_getCarPF(x(3, ii),...
                            lv.v(1, ii),...
                            h,...
                            lv.l,...
                            lv.w,...
                            lv.f.xpos,...
                            lv.f.ypos,...
                            Road.Xg,...
                            Road.Yg);
    else
        lv.f.xpos = nan;
        lv.f.ypos = nan;
        obs.f.pf = 0;
    end
    % exception handling
    if isempty(obs.f.pf)
        obs.f.pf = zeros(size(Road.pf));
    end
    % Combine obstacle pf with road pf to get overall PF
    Road.f.pf = Road.f.pf + obs.f.pf;
    % truncate very large values at road edges
    Road.f.pf(Road.f.pf > 150) = 150;
    Road.f.pf = round(Road.f.pf.*1e4)/1e4;

    sv.reach.pf     = Road.f.pf(sv.reach.min.yi:sv.reach.max.yi, sv.reach.min.xi:sv.reach.max.xi);
    % find minimum pf in given reachable space
    sv.reach.minpf.pf = min(unique(sv.reach.pf));
    % furthest point in reachable space with minimum pf (unique point)
    [sv.reach.minpf.yi, sv.reach.minpf.xi, ~] = find(Road.f.pf(:, 1:sv.reach.max.xi) == sv.reach.minpf.pf, 1, 'last');
    % extract (x,y) coordinate of minpf point
    sv.reach.minpf.y = Road.Yg(sv.reach.minpf.yi, 1);
    sv.reach.minpf.x = Road.Xg(1, sv.reach.minpf.xi);

    % eliminate minor oscillations in lateral target
    sv.reach.minpf.y = max(lw/2,min(3*lw/2,sv.reach.minpf.y));

    %%
    % create reference signal based on minimum PF
    xt(:,ii) = [sv.reach.minpf.y - lw/2;...
                0;...
                sv.reach.minpf.x/h];
        if xt(1,ii) > 1
            xt(1,ii) = lw;
            xt(3,ii) = 0.95*max(vp);
        else
            xt(1,ii) = 0;
        end

    % extract snapshots
    if (rem(tSNAP(ii),5)<=1e-5)
        tempxSnap = [lv.x(1,ii); lw/2; sv.reach.minpf.x + s(1,ii); min(sv.reach.minpf.y,3*lw/2); s(1,ii); tSNAP(ii)];
        xsnap = horzcat(xsnap,tempxSnap);
    end

    % % trajectory planning % %

	% % Start of section: Add time varying collision avoidance constraints % %
        if lv.x(1,ii) - s(:, ii) < 50 && lv.x(1, ii) - s(:, ii) > 50
            % front collision line (ax + by + c = 0)
            fntColl.pt1.x = lv.x(1,ii) - lv.l/2 - (x(3,ii)*h);
            fntColl.pt1.y = lv.y(1,ii);
            fntColl.pt2.x = lv.x(1,ii) - lv.l/2;
            fntColl.pt2.y = lv.y(1,ii) + lv.w/2;
            %
            fntColl.a     = (fntColl.pt1.y - fntColl.pt2.y);
            fntColl.b     = (fntColl.pt2.x - fntColl.pt1.x);
            fntColl.c     = (fntColl.pt1.x*fntColl.pt2.y - fntColl.pt2.x*fntColl.pt1.y);

            % rear collision line (ax + by + c = 0)
            rearColl.pt2.x = lv.x(1,ii) + lv.l/2 + (x(3,ii)*h);
            rearColl.pt2.y = lv.y(1,ii);
            rearColl.pt1.x = lv.x(1,ii) + lv.l/2;
            rearColl.pt1.y = lv.y(1,ii) + lv.w/2;
            %
            rearColl.a     = (rearColl.pt1.y - rearColl.pt2.y);
            rearColl.b     = (rearColl.pt2.x - rearColl.pt1.x);
            rearColl.c     = (rearColl.pt1.x*rearColl.pt2.y - rearColl.pt2.x*rearColl.pt1.y);

            % collision avoidance boundary defined by [ ca.a*x + ca.b*y + ca.c = 0 ]
            if s(:,ii) < lv.x(1,ii) % use front collision constraint
                ca.a = fntColl.a;
                ca.b = fntColl.b;
                ca.c = fntColl.c;
            else % use rear collision constraint
                ca.a = rearColl.a;
                ca.b = rearColl.b;
                ca.c = rearColl.c;
            end
            % constraint VII (collision avoidance consiitraints) - for initial nominal state
            L7 = -[ca.b, 0, 0, zeros(1,Hc*nu), zeros(1,np)];
            c7 = (ca.a*s(:,ii) + ca.c);
            W7 = zeros(1,nx+nx);

            % constraint VIII (collision avoidance constraints) - for predicted states
            L8  = zeros(Hp,nx+(Hc*nu)+np);
            c8  = c7*ones(Hp,1);
            W8  = zeros(Hp,nx+nx);
            for mm = 1:Hp
                xi_pred = zeros(1, nx + Hc*nu);
                for nn = 1:mm
                    xi_pred = xi_pred + [ Phi(nn*nx,:), Gamma(nn*nx,:) ];
                end
                xi_pred = ca.a*ts*xi_pred;
                % replace first row of temp matrix with mm^th predicted lateral position
                eta_pred = ca.b*[ Phi((mm-1)*nx+1,:), Gamma((mm-1)*nx+1,:) ];

                L8(mm,1:nx+(Hc*nu)) = -(xi_pred + eta_pred);
            end
        else
            L7 = []; c7 = []; W7 = [];
            L8 = []; c8 = []; W8 = [];
        end
    % % End of section % %

    % constrained MPC
    tic
    [Xopt, FVAL, exitFlg(ii)] = quadprog( G,...
                                          F*[x(:,ii);xt(:,ii)],...
                                          [],...%[L;L7;L8],...
                                          [],...%[c;c7;c8]+[W;W7;W8]*[x(:,ii);xt(:,ii)],...
                                          [],...
                                          [],...
                                          [],...
                                          [],...
                                          [],...
                                          options );   % constrained MPC
    tictoc(1, ii) = toc;

	% extract optimal input and optimised initial condition of nominal system
    if isempty(Xopt)
        u(:,ii) = zeros(nu,1);

        xbar0   = x(:,ii);
    else
        u(:,ii) = Xopt(nx+1:nx+nu,:);

        xbar0	= Xopt(1:nx,:);
    end
    % control law
    u_rob(:,ii) = u(:,ii) - K_e*( x(1:nx,ii) - xbar0 );

    % generate predicted state trajectories
    xp = Phi*x(:,ii) + Gamma*( Xopt(nx+1:nx+Hc*nu,:) + [-K_e*(x(1:nx,ii)-xbar0);repmat([0;0],Hc-1,1)] );
    traj.y  = zeros(1, Hp);
    traj.vp = zeros(1, Hp);
    traj.x  = zeros(1, Hp);
    for mm = 1:Hp
        traj.y(1, mm)   = xp(1+nx*(mm-1)) + lw/2;
        traj.vp(1, mm)  = xp(3+nx*(mm-1));
        traj.x0         = traj.vp(1, mm)*ts;
        if mm == 1
            traj.x(1, mm) = traj.x0;
        else
            traj.x(1, mm) = traj.x0 + traj.x(1, mm-1);
        end
    end

    % non-linear bicycle model acting as plant
    tsnap   = 0:1e-3:ts;
    acc     = linspace(u_rob(1, ii), u_rob(1, ii), length(tsnap));
    delta   = linspace(u_rob(2, ii), u_rob(2, ii), length(tsnap));
    if ii > 1
        [~, Xp] = ode45(@(t, x) func_nlKinBicycle(t, x, tsnap, acc, delta), tsnap, [s(:, ii); x(:,ii)]);
    else
        [~, Xp] = ode45(@(t, x) func_nlKinBicycle(t, x, tsnap, acc, delta), tsnap, [s(:, ii); x0]);
    end
    % saving vehicle states
    s(:, ii+1) = Xp(end, 1);
    x(:, ii+1) = Xp(end, 2:end)';

    vertices = [bnd.x(bnd.B),bnd.y(bnd.B) + x(1,ii) + lw/2];

    % % UNCOMMENT THE LINES BELOW TO GET PLOT ANIMATIONS % %
    %{
    if lv.x(1,ii) - s(:,ii) < 50 || lv.x(1,ii) - s(:,ii) > -50
        traj.xStack(kk,:) = ((x(3,ii) - lv.v(1,ii))/x(3,ii))*traj.x + s(:,ii) - lv.x(1,ii);
        traj.yStack(kk,:) = traj.y - lw/2;

        traj.xFollow(:,kk) = s(:,ii) - lv.x(1,ii);
        traj.yFollow(:,kk) = x(1,ii);
        traj.yTarget(:,kk) = xt(1,ii) + lw/2;

        kk = kk + 1;
    end
    %}
    %{
    % plot perception and trajectory planning
        figure(10);
            clf;
            subplot(4,1,1);
                contour(Road.Xg, Road.Yg, Road.pf,1e1); hold on;
                    str = sprintf('Surrounding PF at T = %.2f and distance = %.2f',t(ii),lv.x(1,ii)-s(:,ii));
                    title(str);
                plot(sv.block.x, sv.block.y, 'b',...
                     lv.block.x, lv.block.y, 'r');
                plot([0,traj.x(1,1:end)], [x(1,ii+1)+lw/2, traj.y(1,1:end)], 'k', 'LineWidth', 1); hold on;
                plot(Road.Xg(1,sv.reach.minpf.xi),Road.Yg(sv.reach.minpf.yi,1),'m+');
            subplot(4,1,2:4);
                surf(Road.Xg,Road.Yg,Road.pf);
                    xlabel('x co-ordinates [m]');
                    ylabel('y co-ordinates [m]');
                    zlabel('PF(x,y) magnitude [-]');
                    az = -71;
                    el = 38;
                    view(az, el);
    %         subplot(4,1,4);
    %             plot(t(1,1:ii),exitFlg(1,1:ii),'b'); grid on; hold on;
    %                 xlim([0 43.4 ]);
    %                 ylim([-6 2]);
    %}

    %{
    fig = figure(11);
        clf;
        surf(Road.Xg,Road.Yg,Road.pf);
            xlabel('x co-ordinate [m]');
            xlim([rrear rfront]);
            ylabel('y co-ordinate [m]');
            zlabel('PF(x,y) [-]');
            az = -71;
            el = 38;
            view(az, el);
    pause(8e-2);
    %}

    %{
    x0      = 25;
    y0      = 1;
    width   = 18;
    height  = 21;

    figure(12);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[x0,y0,width,height]);
    set(gcf,'PaperPositionMode','auto');
    clf;
        subplot(4,2,1:2);
            contour(Road.Xg, Road.Yg, Road.pf,1e1); hold on;
                str = sprintf('Distance to vehicle = %.2f [m]',lv.x(1,ii)-s(:,ii));
                title(str);
                ylabel('$\eta_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex');
                xlabel('$\xi_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex')
                plot(sv.block.x, sv.block.y, 'b',...
                     lv.block.x, lv.block.y, 'r','LineWidth',2);
            %{
            plot(traj.x, traj.y, 'k--', 'LineWidth', 1); hold on;
            plot(Road.Xg(1,sv.reach.minpf.xi),Road.Yg(sv.reach.minpf.yi,1),'m+', 'LineWidth', 3);
            %}
            patch(vertices(:,1:2),vertices(:,3:4),'y','LineWidth',1,'EdgeColor','y');
                alpha(0.3);
                ylim([-0.2 7.2]);
                xlim([-60 100]);
        subplot(4,2,3:4);
            plot(s(1,1:ii), x(1,1:ii),'b','LineWidth',2);    grid on; hold on;
            plot(lv.x(1,1:ii), lv.y(1,1:ii),'r-.','LineWidth',2);
                legend('SV','LV','Location','NorthEast');
            plot(s(1,ii), x(1,ii),'bd');
            plot(lv.x(1,ii), lv.y(1,ii),'rd');
                xlim([0 1100.5]);
                ylim([-lw/2 nlane*lw]);
                ylabel('$\eta~ \rm{[m]}$','Interpreter','Latex');
                xlabel('$\xi~ \rm{[m]}$','Interpreter','Latex');
        subplot(4,2,5);
            plot(t(1,1:ii),x(3,1:ii),'b'); grid on; hold on;
            plot(t(1,1:ii),lv.v(1,1:ii),'m','LineWidth',1);
            plot(t(1,1:ii),xt(3,1:ii),'g--','LineWidth',1);
            plot(t, ones(1,length(t)).*xmax(3,1), 'r--', 'LineWidth', 1);
            plot(t, ones(1,length(t)).*xmin(3,1), 'r--', 'LineWidth', 1);
                xlim([0 max(t)]);
                ylim([20 35]);
                ylabel('$v ~\rm{[m/s]}$','Interpreter','Latex');
        subplot(4,2,7);
            plot(t(1,1:ii),u_rob(1,1:ii),'b'); grid on; hold on;
            plot(t, ones(1,length(t)).*umax(1,1), 'r--', 'LineWidth', 1);
            plot(t, ones(1,length(t)).*umin(1,1), 'r--', 'LineWidth', 1);
                xlim([0 max(t)]);
                ylim([-2 2]);
                ylabel('$a_{\rm{x}} ~\rm{[m/s^2]}$','Interpreter','Latex');
                xlabel('$\rm{Time ~[s]}$','Interpreter','Latex');
        subplot(4,2,6);
            plot(t(1,1:ii),radtodeg(x(2,1:ii)),'b'); grid on; hold on;
            plot(t, radtodeg(ones(1,length(t)).*xmax(2,1)), 'r--', 'LineWidth', 1);
            plot(t, radtodeg(ones(1,length(t)).*xmin(2,1)), 'r--', 'LineWidth', 1);
                xlim([0 max(t)]);
                ylim([-2.5 2.5]);
                ylabel('$\psi ~\rm{[deg.]}$','Interpreter','Latex');
        subplot(4,2,8);
            plot(t(1,1:ii), radtodeg(u_rob(2,1:ii)),'b'); grid on; hold on;
            plot(t, radtodeg(ones(1,length(t)).*umax(2,1)), 'r--', 'LineWidth', 1);
            plot(t, radtodeg(ones(1,length(t)).*umin(2,1)), 'r--', 'LineWidth', 1);          
                xlim([0 max(t)]);
                ylim([-2 2]);
                ylabel('$\delta_{\rm{f}} ~\rm{[deg.]}$','Interpreter','Latex');
                xlabel('$\rm{Time ~[s]}$','Interpreter','Latex');
    %}

    %{
    x0      = 05;
    y0      = 1;
    width   = 20.8;
    height  = 20.7;

    figure(12);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[x0,y0,width,height/3]);
    set(gcf,'PaperPositionMode','auto');
    clf;
        contour(Road.Xg, Road.Yg, Road.pf,1e1); hold on;
            ylabel('$\eta_{\rm{r}}~ \rm{[m]}$','Interpreter','Latex');
            xlabel('$\xi_{\rm{r}}~ \rm{[m]}$','Interpreter','Latex');
        plot(sv.block.x, sv.block.y, 'b',...
             lv.block.x, lv.block.y, 'r','LineWidth',2);
        plot(Road.Xg(1,sv.reach.minpf.xi),Road.Yg(sv.reach.minpf.yi,1),'m+', 'LineWidth', 3);
    %     plot(cnvxP,'color','y', 'LineWidth', 0.1); alpha(0.2);
        patch(vertices(:,1:2),vertices(:,3:4),'y','LineWidth',1,'EdgeColor','y');
            alpha(0.2);
            grid off
            box on
            ylim([-0.2 7.2]);
            xlim([-60 100]);
    %}
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return
script_plotResults;