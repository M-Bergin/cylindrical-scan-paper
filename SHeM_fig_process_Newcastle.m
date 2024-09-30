function SHeM_fig_process_Newcastle(meas,median_diff,Scalebar_length, sb_pos)
%Script to process images for publication
%MB 22/10/2018

%% Post process the image if needed and plot


%Remove median height difference if needed
if median_diff
    Z=median_diff_removal(meas.image);
else
    Z=meas.image;
end


%figure;imagesc((meas.inputs{10}),(meas.inputs{6}),Z);colormap(gray); axis equal tight
figure;imagesc([meas.start(1), meas.start(1)+meas.image_size(1)],[meas.start(2), meas.start(2)-meas.image_size(2)],Z); colormap gray;axis square equal tight;
set(gca,'YDir','normal');

% if flip==1
%     set(gca,'YDir','normal','XDir','reverse')
% end

%% Add scalebar and label to image
%Shifts on the scalebar from the edge of the image
vert_shift=0.06;
hori_shift=0.1;

y_lims=[meas.start(2), meas.start(2)-meas.image_size(2)];
x_lims=[meas.start(1), meas.start(1)+meas.image_size(1)];


if sb_pos==0
    x_location=x_lims(1)+hori_shift*abs(diff(x_lims));
    y_location=y_lims(2)+vert_shift*abs(diff(y_lims));
elseif sb_pos==1
    x_location=x_lims(2)-hori_shift*abs(diff(x_lims));
    y_location=y_lims(1)-vert_shift*abs(diff(y_lims));
end


hold on

%Add scalebar
if sb_pos==0
    q_h=quiver(x_location,y_location,Scalebar_length,0,'ShowArrowHead','off','AutoScale','off','LineWidth',5,'Color','y');
elseif sb_pos==1
    q_h=quiver(x_location,y_location,-Scalebar_length,0,'ShowArrowHead','off','AutoScale','off','LineWidth',5,'Color','y');
end

% %Add label
% if sb_pos==0
%     t_h=text(x_location+Scalebar_length/2,y_lims(2)+(vert_shift*2.2)*abs(diff(y_lims)),[num2str(Scalebar_length),' ',char(181),'m'],'Color','y','FontSize',24,'HorizontalAlignment','center');
% elseif sb_pos==1
%     t_h=text(x_location-Scalebar_length/2,y_lims(1)-(vert_shift*2)*abs(diff(y_lims)),[num2str(Scalebar_length),' ',char(181),'m'],'Color','y','FontSize',24,'HorizontalAlignment','center');
% end
%% Tidy up the figure

%Clear the axes
axis off


%Clear the greyspace

im_size=size(meas.image);

fig_h=gcf;
ax_h=gca;
%fig_h.PaperPositionMode='auto';

set(ax_h,'units','pixels') % set the axes units to pixels
x = get(ax_h,'position'); % get the position of the axes
set(fig_h,'units','pixels') % set the figure units to pixels
y = get(fig_h,'position'); % get the figure position
set(fig_h,'position',[y(1) y(2) (im_size(2)/im_size(1))*x(4) x(4)])% set the position of the figure to the length and width of the axes
set(ax_h,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels



%Set paper size for printing to pdf

set(fig_h,'Units','Inches');
pos = get(fig_h,'Position');
set(fig_h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])



%% Optionally add box for another image


% %Grab z position if it is in the datafile
% try
%     z_original=meas.initial_pos(3);
% catch
%     z_original=NaN;
% end
% 
% %Load zoomed in image
% load()
% 
% 
% y_lims=meas.inputs{6};
% x_lims=meas.inputs{10};
% try
%     z_zoom=meas.initial_pos(3);
% catch
%     z_zoom=NaN;
% end
% m=0.9887;
% m2=-0.00055;
% 
% 
% %Check difference in z
% z_diff=double(z_original-z_zoom);
% if isnan(z_diff)
% 
% elseif abs(z_diff)>100
%     warning('Different z positions')
%     x_lims=x_lims+m*z_diff;
% end
% 
% hold on
% r_h=rectangle('Position',[x_lims(1) y_lims(1) x_lims(2)-x_lims(1) y_lims(2)-y_lims(1)],'EdgeColor','r','LineWidth',1 );
% 
