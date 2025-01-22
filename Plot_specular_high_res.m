% Plot the high res scan of specular with the straight through scan

%% Create figure

fig_h=figure;
fig_h.Position=[500 218 1120 340];
tiledlayout(1,2);

plt_label_size=25;


%% Plot the straight through direct scan

load('data\An014096.mat')

%Plot the raw data
% nexttile
subplot(1,2,1)

imagesc(angle_sorted/1e6,(det_y_sorted/1e6)-0.279,counts_mat_sorted);
ylabel('\Deltah/mm','FontSize',plt_label_size)
xlabel('\phi/^\circ','FontSize',plt_label_size)
set(gca,'FontSize',plt_label_size,'LineWidth',1)
% title('Direct')
ax1=gca;
% ax=gca;
% ax.Color=[ 0.2422    0.1504    0.6603];
% 
% ylim([0 0.7])

%ylim([-0.1 0.54])
ylim([-0.25 0.25])
set(ax1,'Color',[0.2422    0.1504    0.6603])



%% Plot specular
load('data\An004632.mat')

%Plot the raw data
%nexttile
subplot(1,2,2)
imagesc(angle_sorted/1e6,(det_y_sorted/1e6)-0.14,counts_mat_sorted);
ylabel('\Deltah/mm','FontSize',plt_label_size)
xlabel('\phi/^\circ','FontSize',plt_label_size)

xlim([44.25 51.75])
ylim([-0.25 0.25])
% ylim([-0.35 0.35])
% ylim([-0.16 0.54])
%ylim([-0.25 0.45])
%xticks([45 47 49 51]);
set(gca,'FontSize',plt_label_size,'LineWidth',1)
% title('Specular')
ax2=gca;
temp_position=ax2.Position;



 %% Export result
%exportgraphics(fig_h,'..\specular.pdf','ContentType','vector')

