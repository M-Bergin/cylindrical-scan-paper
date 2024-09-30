function [data, theta_out] = load_an_scans_MB(files_ind,path_name)


fname = [path_name,'/An', num2str(files_ind(1), '%06.f'), '.mat'];
load(fname,'angle_vec')

theta_out = angle_vec'/1e6;



for ind=1:length(files_ind)
    fname = [path_name,'/An', num2str(files_ind(ind), '%06.f'), '.mat'];
    load(fname,'counts')
    
    data(:,ind) = counts;

    
end


%data=data';
end

