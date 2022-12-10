function R = func_getRoadPF(lw,nlane)
    % getRoadPF  Obtain potential field distribution of a given road type
    % R = getRoadPF(lw,nlane)
    %
    % Inputs
    %       lw      : lane width [m]
    %       nlane   : total number of lanes on the road [-]
    %
    % Outputs
    %       R       : structure containing (x,y) coordinates and value of
    %                 potential function for a given road segment at the given points
    %
    % Author   : Shilp Dixit
    % Date     : 04 July 2017
    % Location : 3AB01 UoS

    %%
    rfront = 100;     % sensor range front [m]
    rrear  = -60;     % sensor range rear [m]

    % spatial resolution
    ystp = 0.1;     % lateral position step size [m]
    xstp = 02;      % longitudinal position step size [m]

    % lane potential parameters
    Alane = 36;     % Max height
    sigma = 0.14*lw;% width parameter

    % road potential parameters
    eta   = 3;      % scale

    % lane velocity parameters
    gamma = 0.2;
    vdes  = 27.77;  % desired vehicle speed [m/s]
    vmax  = 36.11;  % maximum vehicle speed [m/s]
    vlane = linspace(vdes,vmax,nlane);  % array of vehicle speeds at each lane [m/s]

    % initialize grid points
    xg = [rrear:xstp:rfront];
    yg = [0:ystp:nlane*lw];
    % create mesh
    [Xg,Yg] = meshgrid(xg,yg);

    % calculate lane division PF
    l.pf = zeros(size(Xg,1),size(Xg,2),nlane-1);
    for ii = 1:nlane-1
        l.pf(:,:,ii) = Alane*exp(-(Yg - ii*lw).^2/(2*sigma^2));
    end

    % calculate road PF
    R.pf = 0*exp(-(Xg.^2+Yg.^2));
    for ii = 0:1
        R.pf = R.pf + (1/2)*eta*1./(Yg - ii*nlane*lw).^2;
    end

    % calculate PF for each lane based on longitudinal velocity
    v.pf = 0*exp(-(Xg.^2+Yg.^2));
    for ii = 1:nlane
        ys = find(Yg(:,1) == (ii-1)*lw,1,'first');  % calculate index of start of lane boundary
        yf = find(Yg(:,1) == ii*lw,1,'first');      % calculate index of end of lane boundary

        v.pf(ys:yf,:) = gamma*(vlane(ii) - vdes);   % assign pf due to speed of the lane
    end

    % calculate PF of road + lane marking
    for ii = 1:nlane-1
        R.pf = R.pf + l.pf(:,:,ii);
    end

    % add lane velocity PF to net PF
    R.pf = R.pf + v.pf;
    R.pf(isinf(R.pf)) = nan;

    R.Xg = Xg;
    R.Yg = Yg;
return