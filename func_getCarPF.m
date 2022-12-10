function O = func_getCarPF(v,vobs,htime,l,w,xpos,ypos,Xg,Yg)
    % getCarPFv4  Obtain potential field distribution of an obstacle
    % [Xg, Yg, pf] = getCarPF(v,vobs,htime,l,w,xpos,ypos,Xg,Yg)
    %
    % Inputs
    %       v       : velocity of subject vehicle [m/s]
    %       vobs    : velocity of obstacle [m/s]
    %       htime   : headway time [s]
    %       l       : lengh of obstacle [m]
    %       w       : width of obstacle [m]
    %       xpos    : x position of obstacle w.r.t source [m]
    %       ypos    : y position of obstacle w.r.t source [m]
    %       Xg      : mesh x-grid points of given space
    %       Yg      : mesh y-grid points of given space
    %
    % Outputs
    %       Xg      : mesh x-grid points of given space
    %       Yg      : mesh y-grid points of given space
    %       pf      : potential function of obstacle
    %
    % Author   : Shilp Dixit
    % Date     : 11 July 2017
    % Location : 3AB01 UoS

    %%
    % spaial discretisation parameters
    ystp = 0.05;    % lateral position step size [m]
    xstp = 02;      % longitudinal position step size [m]

    % car PF parameters
    Acar    = 10;   % Yukawa amplitude
    alpha   = 0.6;  % Yukawa scale
    %%
    % position of obstacle
    xg = [xpos:xstp:xpos+l];
    yg = [ypos - w/2:ystp:ypos + w/2];

    % calculate co-ordinates lying in the virtual polyon
    car.xv = [min(max(Xg(1,:)),max(xg) + vobs*(htime)),...
              max(xg),...
              min(xg),...
              max(min(Xg(1,:)),min(xg) + v*(-htime)),...
              min(xg),...
              max(xg),...
              min(max(Xg(1,:)),max(xg) + vobs*(htime))];
    car.yv = [yg(1,ceil(end/2)),...
              min(yg),...
              min(yg),...
              yg(1,ceil(end/2)),...
              max(yg),...
              max(yg),...
              yg(1,ceil(end/2))];

    % identify points in grid lying inside the polygon (car + 2 virtual triangles)
    [in,~] = inpolygon(Xg,Yg,car.xv,car.yv);
    % converting from logical to double
    in = double(in);
    
    [car.r, car.c] = find(in == 1);
    % grid points occupied by vehicle
    car.xg = Xg(1,car.c);
    car.yg = Yg(car.r,1);

    O.pf = zeros(size(Xg));
    K = zeros(size(Xg,1),size(Xg,2),length(car.xg));
    for ii = 1:length(car.xg)
        K(:,:,ii) = sqrt((car.xg(1,ii)-Xg).^2 + (car.yg(ii,1)-Yg).^2);
    end
    K(K==0) = ystp;

    Kf = min(K,[],3);
    O.pf = (Acar./Kf).*exp(-alpha.*Kf);
    
    O.Xg = Xg;
    O.Yg = Yg;
return;