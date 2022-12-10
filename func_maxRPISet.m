function [h_Fsplot, Fs, h_Fasplot, Fas] = func_maxRPISet( W, Ak, loopCount, epsilon )
    % func_maxRPISet( W, Ak, loopCount, epsilon )
    % used to create the minimum robust positive invariant set for a system
    %
    % [h_Fsplot, Fs, h_Fasplot, Fas] = func_maxRPISet( W, Ak, loopCount, epsilon )
    %
    % Inputs
    %       W               : disturbance set
    %       Ak              : stable system matrix
    %       loopCount       : iteration limit
    %       epsilon         : outer set approximation limit
    %
    % Outputs
    %       h_Fsplot        : plot handle of Fs plot
    %       Fs              : Final set in teration for rpi approximation
    %       h_Fasplot       : plot handle of Fas plot
    %       Fas             : minimal robust positively invariant set for the error
    %
    %	Author   : Shilp Dixit
    %   Date     : 01 December 2017
    %   Location : 17AA03 UoS

    % dummy ouput instead of plot handles
    h_Fsplot = 0;
    h_Fasplot = 0;

    % unit vectors in state/disturbance space
    e = eye( size(Ak,1) );

    for ii = 1:loopCount
        s = ii;

        alpha0 = 0;
        for jj = 1:size(W.A,1)
            alpha0 = max( alpha0,(W.support( Ak^(s)'*W.A(jj,:)' )/W.b(jj)) );
        end
        alphas = min(alpha0,1e4);

        for kk = 1:size(e,2)
            M1 = 0;
            M2 = 0;
            for ll = 1:s
                M1 = W.support( Ak^(ll-1)'*e(:,kk) ) + M1;
                M2 = W.support( -(Ak^(ll-1))'*e(:,kk) ) + M2;
            end
            M = max(M1,M2);
        end
        Ms = max(M);

        if s == 1
            Fs = Ak^(s-1)*W;
        else
            Fs = Ak^(s-1)*W + Fs;
        end

        % loop break condition
        if alphas <= epsilon/(epsilon + Ms)
            %  plot last intermediate set
            %{
            h_Fsplot = plot(Fs, 'color', 'k', 'linewidth', 0.5, 'linestyle', '-.');
            uistack(h_Fsplot,'bottom');
            alpha(0.4);
            %}

            % scaling of set (to obtain approximation of F_infinity RPI set)
            Fas = (1-alphas)^(-1)*Fs;

            % plot final RPI set
            %{
            h_Fasplot = plot(Fas, 'color', 'y', 'linewidth', 0.5);
            uistack(h_Fasplot,'bottom');
            alpha(0.3);
            %}

            break
        elseif ii == loopCount && alphas > epsilon/(epsilon + Ms)
            disp('Set not found. Try changing loopCount and epsilon');

            break
        end
    end

    return