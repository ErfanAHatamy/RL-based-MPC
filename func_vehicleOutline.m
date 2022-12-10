function block = func_vehicleOutline( x_cog, y_cog, yaw, lf, lr, width )
    % func_affineTransform() perform affine tranformation
    %
    % block = func_predictionModel(A,B,Hp,Hc,Q,P,R)
    %
    % Inputs
    %       x_cog   : x-coordinate of CoG of vehicle
    %       y_cog   : y-coordinate of CoG of vehicle
    %       yaw     : heading (yaw) angle of vehicle
    %       lf      : distance CoG to front axle [m]
    %       lr      : distance CoG to rear axle [m]
    %       width   : width of vehicle [m]
    %
    % Outputs
    %       block   : block containing coordinates of extremes of vehicle
    %
    %	Author   : Shilp Dixit
    %   Date     : 06 November 2017
    %   Location : 17AA03 UoS
    
    %%
    % Front left contact patch
    FL.x = (x_cog + lf)*cos(yaw) - (y_cog + width/2)*sin(yaw);
    FL.y = (y_cog + width/2)*cos(yaw) + (x_cog + lf)*sin(yaw);
    % Front right contact patch
    FR.x = (x_cog + lf)*cos(yaw) - (y_cog - width/2)*sin(yaw);
    FR.y = (y_cog - width/2)*cos(yaw) + (x_cog + lf)*sin(yaw);
    % Rear left contact patch
    RL.x = (x_cog - lr)*cos(yaw) - (y_cog + width/2)*sin(yaw);
    RL.y = (y_cog + width/2)*cos(yaw) + (x_cog - lr)*sin(yaw);
    % Rear right contact patch
    RR.x = (x_cog - lr)*cos(yaw) - (y_cog - width/2)*sin(yaw);
    RR.y = (y_cog - width/2)*cos(yaw) + (x_cog - lr)*sin(yaw);
    
    block.x = [FL.x, FR.x, RR.x, RL.x, FL.x];
    block.y = [FL.y, FR.y, RR.y, RL.y, FL.y];