function [Abar, Bbar] = func_nominalModel( modelBank )
    % func_nominalModel( modelBank )
    % obtain the nominal state-space system from a bank of models
    %
    % [Abar, Bbar] = func_nominalModel( modelBank )
    %
    % Inputs
    %       modelBank   : bank of state-space models
    %
    % Outputs
    %       Abar        : system matrix of nominal state-space system
    %       Bbar        : input matrix of nominal state-space system
    %
    %	Author   : Shilp Dixit
    %   Date     : 01 December 2017
    %   Location : 17AA03 UoS

    Asum = zeros( size(modelBank.A{1}) );
    Bsum = zeros( size(modelBank.B{1}) );
    for ii = 1:size(modelBank.A,1)
        Asum = Asum + modelBank.A{ii};
        Bsum = Bsum + modelBank.B{ii};
    end
    Abar = Asum./size(modelBank.A,1);
    Bbar = Bsum./size(modelBank.A,1);

return