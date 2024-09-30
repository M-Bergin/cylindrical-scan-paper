% Script for plotting experimental diffractions scans and comparing them to
% a theoretical prediction
%
% M Bergin
% 28/9/24

%% 

filenames=[4204,4042,4447,3799];

N_files=length(filenames);


%Speed ratio parameters
S=10;
N_lam=50;  %Number of points to calculate wavelength at

filepath=['data\An00',num2str(filenames(1)),'.mat'];
load(filepath)

lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
k_mag=(2*pi)/lambda;

h=det_y_vec_gen/1e6; %in mm
r_arm=6.5; % in mm

% theta=90-asind((h/r)./(sqrt(1+(h/r).^2)));
theta=acosd((h/r_arm)./(sqrt(1+(h/r_arm).^2)));
phi=angle_vec_gen/1e6;

[theta_mat,phi_mat]=meshgrid(theta,phi);

J=(r_arm.^3*cosd(phi_mat))./(h.^2+r_arm.^2).^2;

[h_mat,phi_mat2]=meshgrid(h,phi);
k_X2=k_mag*(r_arm./sqrt(h_mat.^2+r_arm.^2)).*sind(phi_mat2);
k_Y2=k_mag*(h./sqrt(h_mat.^2+r_arm.^2));



I_unsorted_set=NaN*zeros(size(counts_mat_sorted',1),size(counts_mat_sorted',2),N_files);
alpha_unsorted_set=NaN*zeros(N_files,1);

% Load all the data in
for n_files=1:N_files

    filepath=['data\An00',num2str(filenames(n_files)),'.mat'];
    load(filepath)

    % An004528.mat

    I_unsorted_set(:,:,n_files)=counts_mat_sorted';
    alpha_unsorted_set(n_files)=starting_angles(3)/1e6;

end
%%


%Sort the angles
[alpha_sorted_set,ind_sorted]=sort(alpha_unsorted_set,'descend');
I_sorted_set=I_unsorted_set(:,:,ind_sorted);








%%




% Function to approximate the diffraction pattern from a LiF crystal and
% what signal it would produce in a SHeM

%%%%%%%%%%%% Setup of parameters %%%%%%%%%%%

%Parameters
lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
theta_in=45;

alph=0; %Flux in diffuse component





%% Full loading of Boyao data

str=fileread('data/diffrac10001.out');

% Get the start of each row
tkn=regexp(str,'Required number of z grid points');

N_phi=length(tkn);

phi_Boyao_vec=NaN*zeros(N_phi,1);

for n_phi=1:N_phi

    if n_phi==N_phi
        sub_str=str(tkn(n_phi):end);
    else
        sub_str=str(tkn(n_phi):tkn(n_phi+1)-1);
    end

    tkn_n=regexp(sub_str,'n =');
    tkn_n_end=regexp(sub_str(tkn_n:end),'[\n]','once')+tkn_n-2;
    N_rows=str2double(sub_str(tkn_n+3:tkn_n_end));

    tkn_theta=regexp(sub_str,'theta =');
    tkn_theta_end=regexp(sub_str(tkn_theta:end),'[\n]','once')+tkn_theta-2;
    temp=strsplit(sub_str(tkn_theta+6:tkn_theta_end),' ');
    phi_Boyao_vec(n_phi)=str2double(temp{3});


    varnames{n_phi}=matlab.lang.makeValidName(strcat('phi',num2str(phi_Boyao_vec(n_phi))));



    startRow = 7;
    endRow = 8 + N_rows+1;

    % Format for each line of text:
    %   column1: categorical (%C)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%1C%7f%6f%f%[^\n\r]';

    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(sub_str, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    % Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

    % Allocate imported array to column variable names
    VarName1 = dataArray{:, 1};
    b1_full.(varnames{n_phi}) = dataArray{:, 2};
    b2_full.(varnames{n_phi}) = dataArray{:, 3};
    I_full.(varnames{n_phi}) = dataArray{:, 4};


    % Clear temporary variables
    clearvars filename startRow endRow formatSpec fileID dataArray ans;

end


%% Plot the data

fig_h=figure;
fig_h.Position = [100 100 680 435];
t = tiledlayout(2,4,'TileIndexing','columnmajor');
% set(gca,'FontSize',36,'LineWidth',3)
% fig_h.WindowState='fullscreen';
drawnow

% filename='Exp_sim_comparison.gif';
filename='Exp_sim_comparison';
% phi_exp=(45:5:90);
% phi_exp=(45:-5:0);
% phi_exp=(90:-5:45);
% phi_exp=[45,30,15,0];
phi_exp=[0,15,30,45];

h_grid_min=-4;
h_grid_max=4;
h_grid_N=300;

theta_grid_min=-5;
theta_grid_max=62;
theta_grid_N=300;

row_ind_exp=55;
row_ind_sim=2*theta_grid_N/3;

plt_line_colour=[0.8500 0.3250 0.0980];

caption_labels={'(a)','(f)','(b)','(g)','(c)','(h)','(d)','(i)'};
for n_phi=1:N_files

    phi_in=phi_exp(n_phi);%45-59;

    %%%%%%%%%% Creation of scattering distribution in k space %%%%%%%%%%

    k_mag=(2*pi)/lambda;

    %Generate k space grid to plot diffraction pattern on



    % theta_grid_min=-10;
    % theta_grid_max=90;
    % theta_grid_N=200;

    %r_arm=6.5;

    h_grid_vec=linspace(h_grid_min,h_grid_max,h_grid_N);
    theta_grid_vec=linspace(theta_grid_min,theta_grid_max,theta_grid_N);

    [h_grid,theta_grid]=meshgrid(h_grid_vec,theta_grid_vec);

    J_sim=(r_arm.^3*cosd(theta_grid))./(h_grid.^2+r_arm.^2).^2;

    % figure;imagesc(h_grid_vec,theta_grid_vec,theta_grid)
    % figure;imagesc(h_grid_vec,theta_grid_vec,h_grid)

    % k_X=k_mag*r_arm*sind(theta_grid)./sqrt(r_arm^2+h_grid.^2);
    % k_Y=k_mag*h_grid./sqrt(r_arm^2+h_grid.^2);


    % Calculate k_x/k and k_y/k
    k_X_k=r_arm*sind(theta_grid)./sqrt(r_arm^2+h_grid.^2);
    k_Y_k=h_grid./sqrt(r_arm^2+h_grid.^2);

    % figure;imagesc(h_grid_vec,theta_grid_vec,k_X)
    % figure;imagesc(h_grid_vec,theta_grid_vec,k_Y)

    % k_x_max=1.2e11;
    % k_x_N=500;
    %
    % k_y_max=1.2e11;
    % k_y_N=500;
    %
    % k_x_vec=linspace(-k_x_max,k_x_max,k_x_N);
    % k_y_vec=linspace(-k_y_max,k_y_max,k_y_N);
    %
    % [k_X,k_Y]=meshgrid(k_x_vec,k_y_vec);

    %Set k_Z by energy conservation

    k_Z_k=sqrt(1-(k_X_k.^2+k_Y_k.^2));
    imag_inds=imag(k_Z_k)>0;



    %% Import Boyao data

    % Find the data from library

    ind_B=find(phi_Boyao_vec<phi_in+0.5 & phi_Boyao_vec>phi_in-0.5,1);
    b1 = b1_full.(varnames{ind_B});
    b2 = b2_full.(varnames{ind_B});
    I = I_full.(varnames{ind_B});




    %% Create the diffraction pattern
    I_k=zeros(h_grid_N,theta_grid_N);

    %Use speed ratio to generate vector of wavelengths
    lambda_sigma=lambda/(sqrt(2)*S);
    lambda_vec=linspace(lambda-lambda_sigma*3,lambda+3*lambda_sigma,N_lam);

    %% Loop over each wavelength
    for n_lam=1:N_lam
        %
        % %Calculate positions of the diffraction peaks
        % [k_out,G,theta_out,phi_out,N_eff]=diffraction_peak_locations(theta_in,phi_in,lambda_vec(n_lam));

        %Calculate positions of the diffraction peaks
        [k_out,G,theta_out,phi_out,N_eff,N_x,N_y]=diffraction_peak_locations(theta_in,phi_in,lambda_vec(n_lam));

        N_points=size(k_out,1);

        N_y=-N_y;

        %Set the width of the peaks and how they decay
        peak_width=6e8; %

        %Main loop to add in the diffraction pattern by each channel at a time
        for n_point=1:N_points
            %Get intensity
            Boyao_ind= find(b1==N_x(n_point) & b2==N_y(n_point));
            peak_I=I(Boyao_ind);

            if ~isempty(Boyao_ind)
                %Create gaussians with means at values of k_x/k and k_y/k
                I_x_temp=normpdf(k_X_k*((2*pi)/lambda_vec(n_lam)),k_out(n_point,1),peak_width);
                I_y_temp=normpdf(k_Y_k*((2*pi)/lambda_vec(n_lam)),-k_out(n_point,2),peak_width);

                I_k=I_k+normpdf(lambda_vec(n_lam),lambda,lambda_sigma)*peak_I*(I_y_temp.*I_x_temp).*abs(J_sim);

            else
                % disp('Missing peak')
                % disp(N_x(n_point))
                % disp(N_y(n_point))
            end
        end

    end
    %Normalise the diffraction pattern contribution
    I_k=(I_k/(sum(sum(I_k))))*(1-alph);


    %Add in diffuse component
    theta_mat=(atand(sqrt(k_X_k.^2+k_Y_k.^2)./k_Z_k));
    theta_mat(theta_mat~=real(theta_mat))=NaN;

    %Distribution for diffuse scattering
    I_diff=~isnan(theta_mat);
    %Distribution for the solid angle
    %I_diff=1./real(cosd(theta_mat));

    %Normalise the diffuse contribution
    I_diff=(I_diff./(nansum(nansum(I_diff))))*alph;


    %Calculate total diffraction pattern with diffuse component.
    I_tot=I_k+I_diff;


    %Remove the imaginary parts
    I_tot(imag_inds)=0;

    %% Convolve with instrument response function

    r_aperture=0.25;%
    theta_det=2;

    h_spacing=h_grid_vec(2)-h_grid_vec(1);
    theta_spacing=theta_grid_vec(2)-theta_grid_vec(1);

    det_h_vec_1=0:h_spacing:1.5*r_aperture;
    det_h_vec_2=-fliplr(det_h_vec_1);
    det_h_vec=[det_h_vec_2(1:end-1),det_h_vec_1];

    det_theta_vec_1=0:theta_spacing:1.5*theta_det;
    det_theta_vec_2=-fliplr(det_theta_vec_1);
    det_theta_vec=[det_theta_vec_2(1:end-1),det_theta_vec_1];


    [det_H,det_Theta]=meshgrid(det_h_vec,det_theta_vec);


    det_response= ((det_H/r_aperture).^2 + (det_Theta/theta_det).^2) <1;

    %figure;imagesc(det_response)

    I_conv=conv2(I_tot,det_response,"valid");

    % figure;imagesc(theta_grid_vec,h_grid_vec,I_conv'); axis tight
    % ylabel('h/mm')
    % xlabel('\theta/^\circ')
    % set(gca,'FontSize',16,'LineWidth',1)


    %% Move to momentum space

    k_X2_sim=k_mag*(r_arm./sqrt(h_grid.^2+r_arm.^2)).*sind(theta_grid)*1e-10;
    k_Y2_sim=k_mag*(h_grid./sqrt(h_grid.^2+r_arm.^2))*1e-10;

    I_tot_k=I_tot./abs(J_sim);
    I_conv_k=conv2(I_tot_k,det_response,"same");

    if n_phi==1
        I_conv_k_1=I_conv_k;
    end



    %% Plot the result

    set(0,'CurrentFigure',fig_h)
    nexttile
    pcolor(k_X2/1e10,k_Y2/1e10,(I_sorted_set(:,:,n_phi))./(abs(J)))
    ax1=gca;
    shading flat
    axis equal tight

    % xlabel(['k_x/',Ang, '^{-1}']);
    % ylabel(['k_y/',Ang, '^{-1}']);
    % xlabel('k_x /10^{11}m^{-1}')
    % ylabel('k_y /10^{11}m^{-1}');
    set(gca,'FontSize',14,'LineWidth',1)
    xlim_exp=ax1.XLim;
    ylim_exp=ax1.YLim;
    ax1.YDir='reverse';
    set(gca,'xtick',[])
    if n_phi>1
        set(gca,'ytick',[])
    end
    title(['\alpha = ', num2str(phi_in),'^\circ'])
    text(-0.5,-4.5,caption_labels{2*n_phi-1},'Color','white','FontSize',14)

    if n_phi==1
        hold on
        plot(k_X2(:,row_ind_exp)/1e10,k_Y2(:,row_ind_exp)/1e10,'--','LineWidth',1,'Color',plt_line_colour)
    end


    nexttile
    pcolor(k_X2_sim,-k_Y2_sim,I_conv_k)
    ax2=gca;
    shading flat
    axis equal tight
    if n_phi>1
        set(gca,'ytick',[])
    end
    

    % xlh =xlabel(['k_x/',Ang, '^{-1}']);
    % ylh=ylabel(['k_y/',Ang, '^{-1}']);

    % xlh =xlabel('k_x /10^{11}m^{-1}');
    % % xlh.Position=[0.1 1.6 1];
    % ylh=ylabel('k_y /10^{11}m^{-1}');
    ax2.YDir='reverse';

    % %Plot all the peak labels
    % %Calculate positions of the diffraction peaks at average wavelength
    % [k_out,G,theta_out,phi_out,N_eff,N_x,N_y]=diffraction_peak_locations(theta_in,phi_in,lambda);
    % N_points=size(k_out,1);
    % N_1=N_x;
    % N_2=-N_y;
    % for n_points=1:N_points
    %     if min(k_X2_sim(:))<k_out(n_points,1)/1e10 && max(k_X2_sim(:))>k_out(n_points,1)/1e10 && min(k_Y2_sim(:))<-k_out(n_points,2)/1e10 && max(k_Y2_sim(:))>-k_out(n_points,2)/1e10
    % 
    %         if N_1(n_points)>=0
    %             txt_x= num2str(N_1(n_points));
    %         else
    %             txt_x=['\bar{', num2str(-N_1(n_points)),'}'];
    %         end
    % 
    %         if N_2(n_points)>=0
    %             txt_y= num2str(N_2(n_points));
    %         else
    %             txt_y=['\bar{', num2str(-N_2(n_points)),'}'];
    %         end
    % 
    %         txt = ['$(',txt_x,',',txt_y,')$'];
    %         text(k_out(n_points,1)/1e10,k_out(n_points,2)/1e10,txt,'Interpreter','latex')
    %     end
    % end

    set(gca,'FontSize',14,'LineWidth',1)
    % set(gca,'XTick',[-1,-0.5,0,0.5,1]);
    % set(gca,'YTick',[-1,-0.5,0,0.5,1]);
    % xlim([0 0.9])
    % ylim([-0.7 0.7])

    axis equal
    axis tight
    xlim(xlim_exp)
    ylim(ylim_exp)

    text(-0.5,-4.5,caption_labels{2*n_phi},'Color','white','FontSize',14)

    if n_phi==1
        hold on
        plot(k_X2_sim(:,row_ind_sim),k_Y2_sim(:,row_ind_sim),'--','LineWidth',1,'Color',plt_line_colour)
    end
    

    drawnow
    % exportgraphics(gcf,[filename,num2str(n_phi),'.png'],'Resolution',500)

    %     frame = getframe(fig_h);
    %     im = frame2im(frame);
    %         [A,map] = rgb2ind(im,256);
    %
    %         if n_phi == 1
    %             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.15);
    %         else
    %             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.15);
    %         end



end


Ang = char(197);
% xlabel(t,['k_x/',Ang, '^{-1}'],'FontSize',14,'Interpreter','latex')
% ylabel(t,['k_y/',Ang, '^{-1}'],'FontSize',14,'Interpreter','latex')
xlabel(t,'$k_x/ {\mathring{A}}^{-1}$','FontSize',14,'Interpreter','latex')
ylabel(t,'$k_y/{\mathring{A}}^{-1}$','FontSize',14,'Interpreter','latex')
% t.TileSpacing = 'none';
t.TileSpacing = 'tight';
t.Padding = 'tight';

% 







fig_h2=figure;
fig_h2.Position = [780 100 250 435];
t2 = tiledlayout(2,1,'TileIndexing','columnmajor');
nexttile

I_cut_exp=I_sorted_set(:,row_ind_exp,1)./(J(:,row_ind_exp));

plot(k_X2(:,row_ind_exp)/1e10,(I_cut_exp-min(I_cut_exp))./(max(I_cut_exp)-min(I_cut_exp)),'LineWidth',1,'Color',plt_line_colour)
ylabel('$I_{norm}$','Interpreter','latex')
set(gca,'xtick',[])
set(gca,'ytick',[0,1])
set(gca,'FontSize',14,'LineWidth',1)
axis square
title(['\alpha = 0','^\circ'])
text(-0.5,0.9,'(e)','Color','k','FontSize',14)
nexttile


I_first=I_conv_k_1(:,row_ind_sim);
plot(k_X2_sim(:,row_ind_sim),I_first./max(I_first),'LineWidth',1,'Color',plt_line_colour)
% ylabel('$I$','Interpreter','latex')
ylabel('$I_{norm}$','Interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1)
set(gca,'ytick',[0,1])
axis square
% xlabel(t2,['k_x/',Ang, '^{-1}'],'FontSize',14)
xlabel(t2,'$k_x/ {\mathring{A}}^{-1}$','FontSize',14,'Interpreter','latex')

text(-0.5,0.9,'(j)','Color','k','FontSize',14)
ax2=gca;
ax2.XTick=[0,2,4,6,8];


t2.TileSpacing = 'none';
t2.Padding = 'tight';


% exportgraphics(fig_h,'../LiF_sim_comp.pdf','ContentType','vector','Resolution',1000,'BackgroundColor','none')
% exportgraphics(fig_h,'../LiF_sim_comp.svg','ContentType','vector','BackgroundColor','none','Resolution',1000)
% exportgraphics(fig_h2,'../LiF_sim_comp_pt2.eps','ContentType','vector','BackgroundColor','none')

exportgraphics(fig_h,'../LiF_sim_comp.png','Resolution',1000)
exportgraphics(fig_h2,'../LiF_sim_comp_pt2.png','Resolution',1000)
exportgraphics(fig_h2,'../LiF_sim_comp_pt2.eps')



