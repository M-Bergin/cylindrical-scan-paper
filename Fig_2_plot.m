%[1 1 0] orienatation
filename='data\29-Aug-2023_002.txt';

meas=read_MkII_data(filename);

% Rotate to correct orientation
meas.image=rot90(flipud(meas.image),3);

SHeM_fig_process_Newcastle(meas,0,1000, 0)

filepath='data\An003799.mat';
load(filepath)


x(1) = -double(starting_position(1))/1e3;
y(1) = -double(starting_position(2))/1e3;

hold on
plot(x(1),y(1),'.r','MarkerSize',26)

% exportgraphics(gca,'..\29-Aug-2023_002_sb.eps')




lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
k_mag=(2*pi)/lambda;

h=det_y_vec_gen/1e6; %in mm
r=7; % in mm

% theta=90-asind((h/r)./(sqrt(1+(h/r).^2)));
theta=acosd((h/r)./(sqrt(1+(h/r).^2)));
phi=angle_vec_gen/1e6;

[theta_mat,phi_mat]=meshgrid(theta,phi);

J=(r.^3*cosd(phi_mat))./(h.^2+r.^2).^2;

[h_mat,phi_mat2]=meshgrid(h,phi);
k_X2=k_mag*(r./sqrt(h_mat.^2+r.^2)).*sind(phi_mat2);
k_Y2=k_mag*(h./sqrt(h_mat.^2+r.^2));


I_exp=counts_mat_sorted';

fig_h2=figure;
fig_h2.Position=[100 100 650 650];
pcolor(k_X2/1e10,k_Y2/1e10,I_exp./(abs(J)))
ax1=gca;
shading flat
axis equal tight

set(gca,'FontSize',24,'LineWidth',1)
Ang = char(197);
xlabel(['k_x/',Ang, '^{-1}'])
ylabel(['k_y/',Ang, '^{-1}'])
yticks([-4:2:4])

% exportgraphics(gca,'..\An3799.png','Resolution',600)



%Plot the raw data
fig_h=figure;
imagesc(angle_sorted/1e6,det_y_sorted/1e6,counts_mat_sorted);
fig_h.Position = [100 100 400 400];
ylabel('h/mm')
xlabel('\phi/^\circ')
set(gca,'FontSize',24,'LineWidth',1,'YDir','normal')
% exportgraphics(gca,'..\An3799_natural.png','Resolution',600)



