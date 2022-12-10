%% plot results of trajectory planning simulation
% Author   : Shilp Dixit
% Date     : 01 November 2017
% Location : 17AA03 UoS

%% plot results
clear xlim

x0      = 05;
y0      = 1;
width   = 18.8;
height  = 20.7;

figure('Units','centimeters',...
       'Position',[x0,y0,width,height],...
       'PaperPositionMode','auto');
clf;
    subplot(3,2,1:2);
        plot(s(1,:), x(1,:) + lw/2, 'color', [0.6 0.6 0.6]);    hold on;
        plot(lv.x(1,:), lv.y(1,:) + lw/2, 'color', [0 0 0], 'LineStyle', '-.');
        	legend('SV','LV','Location','NorthEast');
            xlabel('$\xi~ \rm{[m]}$','Interpreter','Latex');
            xlim([0 max(s)]);
            ylabel('$\eta~ \rm{[m]}$','Interpreter','Latex');
            ylim([-lw/2 nlane*lw+lw/2]);
        % mark start and end of trajectory
        plot(s(1,1),x(1,1) + lw/2, 's', 'color', [0.6 0.6 0.6],'LineWidth',1.5);
        plot(s(1,end),x(1,end) + lw/2, 'x', 'color', [0.6 0.6 0.6],'LineWidth',1.5);
        plot(lv.x(1,1), lv.y(1,1) + lw/2, 's', 'color', [0 0 0],'LineWidth',1.5);
        plot(lv.x(1,end), lv.y(1,end) + lw/2, 'x', 'color', [0 0 0],'LineWidth',1.5);
        % road and lane marking
        plot(s(1,:),lw*ones(1,length(s(1,:))),'k--','LineWidth',1.5);
        plot(s(1,:),2*lw*ones(1,length(s(1,:))),'k','LineWidth',1.5);
        plot(s(1,:),zeros(1,length(s(1,:))),'k','LineWidth',1.5);
            box off
    %
    ax(1) = subplot(3,2,3);
%         plot(t(1,:),xt(3,:), 'color', [0.1 0.1 0.9], 'LineStyle', '-.', 'LineWidth', 1); hold on; %grid on;
        plot(t(1,:),x(3,:), 'color', [0.6 0.6 0.6]); hold on; %grid on;
        plot(t(1,:),lv.v, 'color', [0 0 0],'LineStyle','-.');
        plot(t, ones(1,length(t)).*xmax(3,1), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
        plot(t, ones(1,length(t)).*xmin(3,1), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
            xlim([0 max(t)]);
            ylim([20 40]);
            ylabel('$v ~\rm{[m/s]}$','Interpreter','Latex');
    ax(2) = subplot(3,2,5);
        plot(t(1,:),u_rob(1,:),'color',[0.6 0.6 0.6]); hold on; %grid on;
        plot(t, ones(1,length(t)).*umax(1,1), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
        plot(t, ones(1,length(t)).*umin(1,1), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
            xlim([0 max(t)]);
            ylim([-2 2]);
            ylabel('$a_{x} ~\rm{[m/s^2]}$','Interpreter','Latex');
            xlabel('$\rm{Time ~[s]}$','Interpreter','Latex');
    ax(3) = subplot(3,2,4);
        plot(t(1,:),radtodeg(x(2,:)),'color',[0.6 0.6 0.6]); hold on; %grid on;
        plot(t, radtodeg(ones(1,length(t)).*xmax(2,1)), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
        plot(t, radtodeg(ones(1,length(t)).*xmin(2,1)), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
            xlim([0 max(t)]);
            ylim([-2.25 2.25]);
            ylabel('$\psi ~\rm{[deg.]}$','Interpreter','Latex');
    ax(4) = subplot(3,2,6);
        plot(t(1,:), radtodeg(u_rob(2,:)),'color',[0.6 0.6 0.6]); hold on; %grid on;
        plot(t, radtodeg(ones(1,length(t)).*umax(2,1)), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
        plot(t, radtodeg(ones(1,length(t)).*umin(2,1)), 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
            xlim([0 max(t)]);
            ylim([-1.5 1.5]);
            ylabel('$\delta_{\rm{f}} ~\rm{[deg.]}$','Interpreter','Latex');
            xlabel('$\rm{Time ~[s]}$','Interpreter','Latex');
    linkaxes(ax,'x');
return
% ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%         export_fig;
% print('simTrajPlan','-depsc2');
% plot state and input set

ssh = figure();
    subplot(1,3,1:2)
        plot(X,'color','k','linewidth', 0.5,...
             Xbar,'color','b','linewidth', 0.5,...
             Xf,'color','r','linewidth', 0.5);
            hold on;
        plot3(x(1,:),x(2,:),x(3,:),'k','LineWidth',2);
            alpha(0.2);
            xlabel('Y [m]');
            ylabel('\psi [rad.]');
            zlabel('v [m/s]');
            h = legend('$\mathcal{X}$','$\bar{\mathcal{X}}$','$\bar{\mathcal{X}}_{t}$');
            set(h,'interpreter','Latex');
    subplot(1,3,3)
        plot(U,'color','c','linewidth', 0.5,...
             Ubar,'color','g','linewidth', 0.5);
            hold on;
        plot(u_rob(1,:),u_rob(2,:),'k','LineWidth',2);
        plot(u(1,:),u(2,:),'r','LineWidth',2);
            xlabel('a_{x} [m/s^2]');
            ylabel('\delta [rad.]');
            h = legend('$\mathcal{U}$','$\bar{\mathcal{U}}$');
            set(h,'interpreter','Latex');
            alpha(0.2);

%%
x0      = 05;
y0      = 1;
width   = 18.8;
height  = 20.7;

figure;
set(gcf,'Units','centimeters');
set(gcf,'Position',[x0,y0,width,2*height/3]);
set(gcf,'PaperPositionMode','auto');
%     subplot(3,1,1:2);
        mesh(Road.Xg,Road.Yg,Road.pf);
            xlabel('$\xi_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex');
            xlim([rrear rfront]);
            ylabel('$\eta_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex');
            zlabel('$U_{\rm{r}}(w_{\rm{r}}) ~[-]$','Interpreter','Latex');
            az = -71;
            el = 38;
            view(az, el);
    subplot(3,1,3);
        contour(Road.Xg, Road.Yg, Road.pf,1e1); hold on;
            str = sprintf('Distance to LV = %.2f [m]',lv.x(1,ii)-s(:,ii));
            title(str);
            xlabel('$\xi_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex');
            ylabel('$\eta_{\rm{r}} ~\rm{[m]}$','Interpreter','Latex');
            plot(sv.block.x, sv.block.y, 'b',...
                 lv.block.x, lv.block.y, 'r','LineWidth',2);
        plot(traj.x, traj.y, 'k--', 'LineWidth', 1); hold on;
        plot(Road.Xg(1,sv.reach.minpf.xi),Road.Yg(sv.reach.minpf.yi,1),'m+', 'LineWidth', 3);
        plot(cnvxP,'color','y', 'LineWidth', 0.1); alpha(0.2);
            ylim([0 7]);
            xlim([-60 100]);
return

%{
figure('Units','centimeters',...
       'Position',[x0,y0,1.8*width,height/3],...
       'PaperPositionMode','auto'); box off;
clf;
    p1 = plot(traj.xFollow, traj.yFollow, 'Color', [0 0 1], 'LineStyle', '-', 'LineWidth', 2); hold on;
    p2 = plot([-lv.l/2 -lv.l/2 lv.l/2 lv.l/2 -lv.l/2], lv.block.y - lw/2, 'r','LineWidth',3);
for kk = 1:length(traj.xFollow)
    p3 = plot(traj.xStack(kk,:), traj.yStack(kk,:), 'Color', [0.2 0.8 0.2], 'LineWidth', 1); hold on;
end
    xlim([-100 90]);
    ylim([-lw 2.5*lw]);
    ylabel('$Y_{\textrm{o}} ~[m]$','Interpreter','Latex');
    xlabel('$X_{\textrm{o}} ~[m]$','Interpreter','Latex');
    legend([p2 p1 p3], 'LV', 'SV actual trajectory','planned trajectories','Location','NorthWest');
    plot(traj.xFollow,(1/2)*lw*ones(1,length(traj.xFollow)),'k--','LineWidth',1);
    plot(traj.xFollow,(3/2)*lw*ones(1,length(traj.xFollow)),'k','LineWidth',1);
    plot(traj.xFollow,(-1/2)*lw*ones(1,length(traj.xFollow)),'k','LineWidth',1);
    box off
%}
%     export_fig;
%     plot(traj.xFollow,traj.yTarget,'c','LineWidth',1);
% print('simOTManv','-depsc2');

%%
% close all
P = Polyhedron([0,0; 53.12,4.086; 53.33,0; 53.12,-4.086]);

figure('Units','centimeters',...
       'Position',[x0,y0,1.8*width,height/4],...
       'PaperPositionMode','auto');
        plot(P,'color','c'); alpha(0.5); hold on;
        plot(sv.block.x,sv.block.y-lw/2,'b','LineWidth',2);
            ylabel('$\eta_{\textrm{v}}~ \rm{[m]}$','Interpreter','Latex');
            xlabel('$\xi_{\textrm{v}}~ \rm{[m]}$','Interpreter','Latex');
% print('reachSet','-dpng');

figure('Units','centimeters',...
       'Position',[x0,y0,1.8*width,height/3],...
       'PaperPositionMode','auto');
        hold on;
        cc = hsv(length(xsnap(1,:)));
        for ii = 1:length(xsnap(1,:))
            h(ii) = plot(xsnap(5,ii),xsnap(2,ii),'s','color',cc(ii,:),'LineWidth',1.5,'MarkerSize',09);
            plot(xsnap(3,ii),xsnap(4,ii),'x','color',cc(ii,:),'LineWidth',1.5,'MarkerSize',09);
            plot(xsnap(1,ii),xsnap(2,ii),'d','color',cc(ii,:),'LineWidth',1.5,'MarkerSize',09);
        end
            ylabel('$\eta_{\textrm{v}}~ \rm{[m]}$','Interpreter','Latex');
            xlabel('$\xi_{\textrm{v}}~ \rm{[m]}$','Interpreter','Latex');
            xlim([-5 1200]);
            ylim([-1 2.5*lw]);
            labels = {'$t=0~s$',...
                      '$t=5~s$',...
                      '$t=10~s$',...
                      '$t=15~s$',...
                      '$t=20~s$',...
                      '$t=25~s$',...
                      '$t=30~s$',...
                      '$t=35~s$',...
                      '$t=40~s$'};
        legend(h, labels(1:length(h)),'Interpreter','Latex','Location','NorthEastOutside');
        plot(xsnap(5,:),lw*ones(1,length(xsnap(1,:))),'k--','LineWidth',1.5);
        plot(xsnap(5,:),2*lw*ones(1,length(xsnap(1,:))),'k','LineWidth',1.5);
        plot(xsnap(5,:),zeros(1,length(xsnap(1,:))),'k','LineWidth',1.5);
%     export_fig;
% print('otSnapshots','-dpng');