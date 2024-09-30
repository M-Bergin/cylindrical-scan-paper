
% Load 2D-Cylindrical scans

% files_ind = [3541, 3550, 3559, 3568, 3577, 3586, 3595, 3604, 3613, 3622, 3631, 3640, 3649, 3658, 3667, 3676, 3685, 3694];
%files_ind = [2787, 2778, 2796, 2769, 2760, 2751, 2742, 2733, 2724, 2715, 2706];

files_ind = 4684:4701;
[data, thetas] = load_an_scans_MB(files_ind, './data');



%Temp = [32, 40, 50, 60, 80, 100, 120, 140, 160, 180, 200];
Temp = [200, 180, 160, 140, 120, 100, 80, 60, 40, 32, 50, 70, 90, 110, 130, 150, 170, 190];



Maximum = max(data,[],1);
Minimum = min(data,[],1);

Maximum = reshape(Maximum,[],1);
Minimum = reshape(Minimum,[],1);

Peak_signal = Maximum - Minimum;
frac_signal = Peak_signal/Peak_signal(10);
Ln_signal = log(frac_signal);
Temp_kelvin = Temp + 273.15;



%%


% fig_h=figure;
% 
% filename = '3D Sample Temperature.gif'; % Specify the output file name
% N=length(files_ind);
% 
% 
% for n=1:N
%     
%     surf(data(:,:,n))
%     zlim([0 15000])
%     title(['Sample Temperature = ',sprintf('%.1f',Temp(n)),'^{\circ}C'])
%     set(gca,'FontSize',12,'LineWidth',1)
%     
% %     frame = getframe(fig_h);
% %     im = frame2im(frame);
% %     [A,map] = rgb2ind(im,256);
% %     
% %     if n == 1
% %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% %     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% %     end
%     
%     pause(0.5)
% 
% end





%%

Middle = data(5,:,:);
N_theta=length(thetas);

% fig_j=figure;
% ax1=gca;
% plot(ax1,NaN,NaN,'.','LineWidth',1,'MarkerSize',14)

filename = '2D Sample Temperature.gif'; % Specify the output file name
N=length(files_ind);

a_vec=NaN*zeros(N,1);
a_vec_error=NaN*zeros(N,1);
mu_vec=NaN*zeros(N,1);
mu_vec_error=NaN*zeros(N,1);
d_vec=NaN*zeros(N,1);
d_vec_error=NaN*zeros(N,1);
fit_list=cell(N,1);
for n=1:N
    
    data_temp=data(:,n);%reshape(Middle(:,:,n),[N_theta,1]);


    % Fit: 'untitled fit 1'.
    [xData, yData] = prepareCurveData( thetas, data_temp );

    % Set up fittype and options.
    ft = fittype( 'a*(normcdf(x,mu,sigma)-normcdf(x,mu+d,sigma))+m*x+c', 'independent', 'x', 'dependent', 'y' );
    excludedPoints = (xData < 41) | (xData > 55);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 -Inf 0 -Inf 42 0];
    opts.StartPoint = [1000 -2000 4 100 44 0.01];
    opts.Upper = [Inf Inf 5 Inf 50 Inf];
    opts.Exclude = excludedPoints;

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    
    %Save outputs
    error=confint(fitresult,0.6827);
    
    a_vec(n)=fitresult.a;
    a_vec_error(n)=diff(error(:,1))/2;

    mu_vec(n)=fitresult.mu;
    mu_vec_error(n)=diff(error(:,5))/2;

    d_vec(n)=fitresult.d;
    d_vec_error(n)=diff(error(:,3))/2;

    fit_list{n}=fitresult;



    fitted_data=fitresult(thetas);
    fig_j=figure;
    ax1=gca;
    plot(ax1,thetas,data(:,n),'.','LineWidth',1,'MarkerSize',14)
    %ylim([0 15000])
    hold on
    plot(ax1,thetas,fitted_data,'LineWidth',1)
    title(['Sample Temperature = ',sprintf('%.1f',Temp(n)),'^{\circ}C'])
    ylabel('Counts')
    xlabel('\theta/^{\circ}')
    set(gca,'FontSize',12,'LineWidth',1)
    
%     frame = getframe(fig_j);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
%     
%     if n == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     end
    
    %pause(1)

end



%% Fit straight line to data with weighted fit
ln_a_vec=log(a_vec);
ln_a_vec_error=(a_vec_error./a_vec);
weights=1./ln_a_vec_error.^2;



N_points=length(a_vec);

x_bar=sum(Temp_kelvin'.*weights)/sum(weights);
y_bar=sum(ln_a_vec.*weights)/sum(weights);
x2_bar=sum(Temp_kelvin'.^2.*weights)/sum(weights);
xy_bar=sum(ln_a_vec.*Temp_kelvin'.*weights)/sum(weights);


m_hat=(xy_bar-(x_bar.*y_bar))/(x2_bar-x_bar.^2);
c_hat=y_bar-m_hat*x_bar;


residuals=ln_a_vec-(m_hat*Temp_kelvin'+c_hat);

sigma2=(1/(N_points-2))*sum(weights.*residuals.^2)/sum(weights);

sigma_m2=sigma2/(N_points*(x2_bar-x_bar.^2));


fprintf('DWF = (%0.4g +- %0.2g) x 10^(-3)\n',-m_hat*1000,sqrt(sigma_m2)*1000);


%% Plot the result
T_plot=linspace(min(Temp_kelvin),max(Temp_kelvin),100);
I_fitted=m_hat*T_plot+c_hat;


figure;errorbar(Temp_kelvin,ln_a_vec,ln_a_vec_error,'.','MarkerSize',12,'Linewidth',1)
hold on
plot(T_plot,I_fitted,'LineWidth',1.5)
xlabel('T/K')
ylabel('$ln(I)$','Interpreter','latex')
% ylabel('ln(I)')
% xlim([275,525])
xlim([300,480])
ylim([6.5 7.5])
yticks([6.6:0.2:7.4])
set(gca,'Linewidth',1,'FontSize',12)


% plot(Temp_kelvin(n_plot),ln_a_vec(n_plot),'x')
annotation('arrow',[0.38,0.43],[0.5,0.58])


% create smaller axes in top right, and plot on it
% ax2=axes('Position',[.6 .6 .3 .3]);
ax2=axes('Position',[.15 .15 .3 .3]);
box on

% Take one of the data points to highlight
n_plot=6;
data_temp=data(:,n_plot);
fitted_data=fit_list{n_plot}(thetas);

plt_colour=[0 0.4470 0.7410];

plot(ax2,thetas(~excludedPoints),data(~excludedPoints,n_plot),'.','LineWidth',1,'MarkerSize',8,'Color',plt_colour)
hold on
plot(ax2,thetas(~excludedPoints),fitted_data(~excludedPoints),'LineWidth',1.5)
plot(ax2,thetas(excludedPoints),data(excludedPoints,n_plot),'.','LineWidth',1,'MarkerSize',8,'Color',[127, 158, 179]/255)


% ylabel(ax2,'Counts / 10^3')
ylabel(ax2,'Counts')
xlabel(ax2,'\theta/^{\circ}')
set(ax2,'FontSize',12,'LineWidth',1)
ax2.YAxisLocation='right';
ax2.XAxisLocation='top';
% ax2.YAxis.Exponent = 3;


% exportgraphics(gcf,'../DWF_fig_v2.pdf')






%% Can instead fit directly to the data
% weights=1./a_vec_error.^2;

% 
% % Fit: 'untitled fit 1'.
% [xData, yData, weights_1] = prepareCurveData( Temp_kelvin, a_vec, weights );
% 
% % Set up fittype and options.
% ft = fittype( 'exp1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [3000 -0.001];
% opts.Weights = weights_1;
% 
% % Fit model to data.
% [fitresult_DWF, gof] = fit( xData, yData, ft, opts );
% 
% 
% 
% %Plot the final result 
% 
% T_plot=linspace(min(Temp_kelvin),max(Temp_kelvin),100);
% I_fitted=fitresult_DWF(T_plot);
% 
% figure;errorbar(Temp_kelvin,a_vec,a_vec_error,'.','MarkerSize',12,'Linewidth',1)
% hold on
% plot(T_plot,I_fitted,'LineWidth',1)
% xlabel('T/K')
% ylabel('Specular counts')
% 
% set(gca,'Linewidth',1,'FontSize',12)
% 
% error_DWF=confint(fitresult_DWF,0.6827);
% 
% DWF=fitresult_DWF.b;
% DWF_error=diff(error_DWF(:,2))/2;
% 
% fprintf('DWF = (%0.4g +- %0.2g) x 10^(-3)\n',-DWF*1000,DWF_error*1000);
% 
% % figure;errorbar(Temp_kelvin, d_vec,d_vec_error,'x')
% % figure;errorbar(Temp_kelvin, mu_vec,mu_vec_error,'x')
% 
% 
% % 
