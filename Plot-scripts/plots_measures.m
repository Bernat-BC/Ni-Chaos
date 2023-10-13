% This code can be used after running measures.m to plot the densities of
% interest. 

close all

clf(figure(1));                                 
figure(1);
plot(Xs, qem/trapz(Xs, qem), 'k-',...
     Xs, M_vec_exp, 'k.', 'MarkerSize', 35);
ylim([0 max(qem/trapz(Xs, qem))]);
% legend('$q$', 'Interpreter', 'latex', 'Fontsize', 30);
xlim([0.1 1]);
xticks([0.25 0.5 0.75 1]);
set(gca,'XTickLabel',[0.25 0.5 0.75 1], 'FontSize', 30);
title('quasi-ergodic density, $q$ ($\varepsilon =$ ' + string(eps) + ')',...
            'Interpreter', 'latex', 'FontSize', 30);

clf(figure(2));                                 
figure(2);
plot(Xs, qsm/trapz(Xs, qsm), 'k-');
ylim([0 max(qsm/trapz(Xs, qsm))]);
% legend('$m$', 'Interpreter', 'latex', 'Fontsize', 30);
xlim([0.1 1]);
xticks([0.25 0.5 0.75 1]);
set(gca,'XTickLabel',[0.25 0.5 0.75 1], 'FontSize', 30);
title('quasi-stationary density, $m$ ($\varepsilon =$ ' + string(eps) + ')', 'Interpreter', 'latex', 'FontSize', 30);
 
clf(figure(3));                                 
figure(3);
plot(Xs, ev_cond_right/trapz(Xs, ev_cond_right), 'k-');
ylim([0 max(ev_cond_right/trapz(Xs, ev_cond_right))]);
% legend('$v$', 'Interpreter', 'latex', 'Fontsize', 30);
xlim([0.1 1]);
xticks([0.25 0.5 0.75 1]);
set(gca,'XTickLabel',[0.25 0.5 0.75 1], 'FontSize', 30);
title('right-eigenfunction, $v$ ($\varepsilon =$ ' + string(eps) + ')',...
                'Interpreter', 'latex', 'FontSize', 30);


% Uncomment if you want to plot the statinoary density

% clf(figure(4));
% figure(4);
% plot(Xs, sm/trapz(Xs, sm), 'b-');
% title('stationary density, $p$ ($\varepsilon =$ ' + string(eps) + ')',...
%                 'Interpreter', 'latex', 'Fontsize', 25);
% xlim([0.1 1]);
% xticks([0.25 0.5 0.75 1]);
% set(gca,'XTickLabel',[0.25 0.5 0.75 1], 'FontSize', 30);
% legend('p', 'Fontsize', 15);


% Uncomment if you want to plot all functions/densities together

% clf(figure(5));                                                          
% figure(5);
% plot(Xs, qem/trapz(Xs, qem), 'b-',...
%      Xs, qsm/trapz(Xs, qsm), 'r-',...
%      Xs, ev_cond_right/trapz(Xs, ev_cond_right), 'k-',...
%      Xs, M_vec_exp, 'k.',...
%      Xs, M_vec_con, 'g.', 'MarkerSize', 10);
% title('Densities on Expanding region ($\varepsilon =$ ' +...
%          string(eps) + ')', 'Interpreter', 'latex','Fontsize', 25);
% xlim([0.1 1]);
% xticks([0.25 0.5 0.75 1]);
% set(gca,'XTickLabel',[0.25 0.5 0.75 1], 'FontSize', 30);
% legend('QEM', 'QSM', '$\nu(x)$', '$M_\mathrm{exp}$',...
%           '$M_\mathrm{cont}$', 'Interpreter', 'latex', 'Fontsize', 15);