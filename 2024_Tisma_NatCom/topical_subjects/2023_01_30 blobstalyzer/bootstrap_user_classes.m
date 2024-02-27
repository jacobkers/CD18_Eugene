function [hist_counts, h_Nboot, hist_perc, hPERC_err]=bootstrap_user_classes(user_data,binax, Nbootsample)
%user data may look like: count_type 1, count_type 2 .... count_type n
%for processing, labels are just 1--n
%per type, there are N_u counts with N_u the number of users.

%% Histogram settings and bootstrapping
N_classes=length(user_data);
hist_counts = hist(user_data,binax);
%resample the data:
[~,bootsam_idx] = bootstrp(Nbootsample,@mean,user_data);
%build the bootstrap-repeated data sets:
boot_samples=user_data(bootsam_idx);
%build repeated histograms and assess the eror per bin:
h_Nboot = hist(boot_samples,binax);
hist_perc = 100*hist_counts/(sum(hist_counts));
hPERC_err=2* 100*std(h_Nboot')/sum(hist_counts);  %95%