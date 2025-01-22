% Script for plotting experimental diffractions scans and comparing them to
% a theoretical prediction
%
% M Bergin
% 10/1/25

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

    I_unsorted_set(:,:,n_files)=counts_mat_sorted'-min(counts_mat_sorted(:));
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
n_phi=4;
phi_in=phi_exp(n_phi);

if n_phi==4
    mid_ind_exp=44;
elseif n_phi==1
    mid_ind_exp=41;
end

%% Determine the FWHM of the two peaks in experiment

k_mid_exp=k_X2(:,mid_ind_exp);
I_mid_exp=(I_sorted_set(:,mid_ind_exp,n_phi))./(abs(J(:,mid_ind_exp)));
I_mid_norm_exp=I_mid_exp/max(I_mid_exp);

scale_factor=0.5;

% Get the data indices
diff_1_inds_exp=find(k_mid_exp>3e10 & k_mid_exp<7e10);
diff_2_inds_exp=find(k_mid_exp>0 & k_mid_exp<3e10);

% Get the actual data
k_1_exp=k_mid_exp(diff_1_inds_exp);
I_1_exp=I_mid_exp(diff_1_inds_exp);

k_2_exp=k_mid_exp(diff_2_inds_exp);
I_2_exp=I_mid_exp(diff_2_inds_exp);

% Get FWHM of peak 1
max_I_1_exp=max(I_1_exp);
min_I_exp=min(I_mid_exp);

halfmax_I_1_exp=(max(I_1_exp)-min(I_mid_exp))*scale_factor+min(I_mid_exp);

leftIndex_exp = find(I_1_exp >= halfmax_I_1_exp, 1, 'first');
rightIndex_exp = find(I_1_exp >= halfmax_I_1_exp, 1, 'last');

left_k1_exp=k_1_exp(leftIndex_exp-1)+((halfmax_I_1_exp-I_1_exp(leftIndex_exp-1))/(I_1_exp(leftIndex_exp)-I_1_exp(leftIndex_exp-1)))*(k_1_exp(leftIndex_exp)-k_1_exp(leftIndex_exp-1));
right_k1_exp=k_1_exp(rightIndex_exp)+((halfmax_I_1_exp-I_1_exp(rightIndex_exp))/(I_1_exp(rightIndex_exp+1)-I_1_exp(rightIndex_exp)))*(k_1_exp(rightIndex_exp+1)-k_1_exp(rightIndex_exp));

fwhm_1_exp = right_k1_exp - left_k1_exp;

% Get FWHM of peak 2
min_I2_exp=min(I_mid_exp((k_mid_exp>1e10 & k_mid_exp<4e10)));
halfmax_I_2_exp=(max(I_2_exp)-min_I2_exp)*scale_factor+min_I2_exp;

leftIndex2_exp = find(I_2_exp >= halfmax_I_2_exp, 1, 'first');
rightIndex2_exp = find(I_2_exp >= halfmax_I_2_exp, 1, 'last');

left_k2_exp=k_2_exp(leftIndex2_exp-1)+((halfmax_I_2_exp-I_2_exp(leftIndex2_exp-1))/(I_2_exp(leftIndex2_exp)-I_2_exp(leftIndex2_exp-1)))*(k_2_exp(leftIndex2_exp)-k_2_exp(leftIndex2_exp-1));
right_k2_exp=k_2_exp(rightIndex2_exp)+((halfmax_I_2_exp-I_2_exp(rightIndex2_exp))/(I_2_exp(rightIndex2_exp+1)-I_2_exp(rightIndex2_exp)))*(k_2_exp(rightIndex2_exp+1)-k_2_exp(rightIndex2_exp));

fwhm_2_exp = right_k2_exp - left_k2_exp;


% figure; plot(k_mid_exp,(I_sorted_set(:,mid_ind_exp,n_phi))./(abs(J(:,mid_ind_exp))))
% xline(left_k1_exp)
% xline(right_k1_exp)
% yline(halfmax_I_1_exp)
% xline(left_k2_exp)
% xline(right_k2_exp)
% xline(k2_mid_exp)
% yline(halfmax_I_2_exp)

%%%%%%%%%% Creation of scattering distribution in k space %%%%%%%%%%

k_mag=(2*pi)/lambda;

%Generate k space grid to plot diffraction pattern on

h_grid_vec=linspace(h_grid_min,h_grid_max,h_grid_N);
theta_grid_vec=linspace(theta_grid_min,theta_grid_max,theta_grid_N);

[h_grid,theta_grid]=meshgrid(h_grid_vec,theta_grid_vec);

J_sim=(r_arm.^3*cosd(theta_grid))./(h_grid.^2+r_arm.^2).^2;

% Calculate k_x/k and k_y/k
k_X_k=r_arm*sind(theta_grid)./sqrt(r_arm^2+h_grid.^2);
k_Y_k=h_grid./sqrt(r_arm^2+h_grid.^2);

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

% Loop over different speed ratios
N_S=43;
S_vec=linspace(4,25,N_S);

fwhm_1=NaN*zeros(N_S,1);
fwhm_2=NaN*zeros(N_S,1);

for n_S=1:N_S

    I_k=zeros(h_grid_N,theta_grid_N);

    %Use speed ratio to generate vector of wavelengths
    lambda_sigma=lambda/(sqrt(2)*S_vec(n_S));
    lambda_vec=linspace(lambda-lambda_sigma*3,lambda+3*lambda_sigma,N_lam);

    %% Loop over each wavelength
    for n_lam=1:N_lam

        %Calculate positions of the diffraction peaks
        [k_out,G,theta_out,phi_out,N_eff,N_x,N_y]=diffraction_peak_locations(theta_in,phi_in,lambda_vec(n_lam));

        N_points=size(k_out,1);

        N_y=-N_y;

        %Set the width of the peaks and how they decay
        peak_width=5e8; %

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


    %% Move to momentum space

    k_X2_sim=k_mag*(r_arm./sqrt(h_grid.^2+r_arm.^2)).*sind(theta_grid)*1e-10;
    k_Y2_sim=k_mag*(h_grid./sqrt(h_grid.^2+r_arm.^2))*1e-10;

    I_tot_k=I_tot./abs(J_sim);
    I_conv_k=conv2(I_tot_k,det_response,"same");

    if n_phi==1
        I_conv_k_1=I_conv_k;
    end

    I_mid=I_conv_k(:,length(I_conv_k)/2);
    I_mid_norm = I_mid./max(I_mid);
    k_mid=k_X2_sim(:,length(I_conv_k)/2);

    %% Determine the FWHM of the two peaks

    % Get the data indices
    diff_1_inds=find(k_mid>3.3 & k_mid<7);
    diff_2_inds=find(k_mid>-0.5 & k_mid<3.3);

    % Get the actual data
    k_1=k_mid(diff_1_inds);
    I_1=I_mid_norm(diff_1_inds);

    k_2=k_mid(diff_2_inds);
    I_2=I_mid_norm(diff_2_inds);

    % Get FWHM of peak 1
    halfmax_I_1=max(I_1)*scale_factor;

    leftIndex = find(I_1 >= halfmax_I_1, 1, 'first');
    rightIndex = find(I_1 >= halfmax_I_1, 1, 'last');

    left_k1=k_1(leftIndex-1)+((halfmax_I_1-I_1(leftIndex-1))/(I_1(leftIndex)-I_1(leftIndex-1)))*(k_1(leftIndex)-k_1(leftIndex-1));
    right_k1=k_1(rightIndex)+((halfmax_I_1-I_1(rightIndex))/(I_1(rightIndex+1)-I_1(rightIndex)))*(k_1(rightIndex+1)-k_1(rightIndex));

    fwhm_1(n_S) = right_k1 - left_k1;

    % Get FWHM of peak 2
    halfmax_I_2=max(I_2)*scale_factor;

    leftIndex2 = find(I_2 >= halfmax_I_2, 1, 'first');
    rightIndex2 = find(I_2 >= halfmax_I_2, 1, 'last');

    left_k2=k_2(leftIndex2-1)+((halfmax_I_2-I_2(leftIndex2-1))/(I_2(leftIndex2)-I_2(leftIndex2-1)))*(k_2(leftIndex2)-k_2(leftIndex2-1));
    right_k2=k_2(rightIndex2)+((halfmax_I_2-I_2(rightIndex2))/(I_2(rightIndex2+1)-I_2(rightIndex2)))*(k_2(rightIndex2+1)-k_2(rightIndex2));

    fwhm_2(n_S) = right_k2 - left_k2;

    % figure;plot(k_mid_exp/1e10,I_mid_norm_exp)
    % hold on
    % plot(k_mid, I_mid_norm)
    % xline(left_k1)
    % xline(right_k1)
    % yline(halfmax_I_1)
    % xline(left_k2)
    % xline(right_k2)
    % yline(halfmax_I_2)
    % title(S_vec(n_S))
    %
    % xline(left_k1_exp/1e10)
    % xline(right_k1_exp/1e10)
    % yline(halfmax_I_1_exp/max(I_mid_exp))
    % xline(left_k2_exp/1e10)
    % xline(right_k2_exp/1e10)
    % xline(k2_mid_exp/1e10)
    % yline(halfmax_I_2_exp/max(I_mid_exp))
    % drawnow



end
%% Plot the result


fig_h_S=figure;
hold on
plot(S_vec,fwhm_2,'LineWidth',1)
xlabel('S')
Ang = char(197);
ylabel(['FWHM/',Ang, '^{-1}'])


% Add lines for experimental results

SIndex2 = find(fwhm_2 <= fwhm_2_exp/1e10, 1, 'first');
S2=S_vec(SIndex2-1)+((fwhm_2_exp/1e10-fwhm_2(SIndex2-1))/(fwhm_2(SIndex2)-fwhm_2(SIndex2-1)))*(S_vec(SIndex2)-S_vec(SIndex2-1));

ax1=gca;

plot([ax1.XLim(1),S2,S2],[fwhm_2_exp/1e10, fwhm_2_exp/1e10, ax1.YLim(1)],'k--','LineWidth',1)

set(gca,'FontSize',14,'LineWidth',1)


% exportgraphics(fig_h_S,'../S.eps')



