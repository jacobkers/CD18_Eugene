function A030_Bananalyzer_Post_Process(init)
 if ~strcmp(init.exp_label, 'BSG217_oufti')
     load([init.datapath_out, init.exp_label, '_shape_data.mat']);
 else  %correct old naming
     load([init.datapath_out, 'results_raw_clicked_ori.mat']);
     shape_data=wrapup_raw;
end

%% Statistics:
%set limits for categorization
lim_circ=0.95;
lim_symm=0.95;

%wrapup: contour_circ contour_symm
circ=shape_data(:,1);
symm=shape_data(:,2);
usr=shape_data(:,3);

%judge the categories:
is_bananas=find(usr== 1);
is_compact=find(usr== 2);
is_pancake=find(usr== 3);
is_multilobe=find(usr== 4);

%plot that is solely based on user judgment:
figure(174);
%scatter
plot(circ(is_bananas),symm(is_bananas),'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'y'); hold on;
plot(circ(is_compact),symm(is_compact),'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
plot(circ(is_pancake),symm(is_pancake),'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
plot(circ(is_multilobe),symm(is_multilobe),'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'm');
%mean
figure(1);
plot(median(circ(is_bananas)),    median(symm(is_bananas)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'y'); 
plot(median(circ(is_compact)),    median(symm(is_compact)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot(median(circ(is_pancake)),    median(symm(is_pancake)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
plot(median(circ(is_multilobe)),   median(symm(is_multilobe)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
axis equal

figure(174);
plot(median(circ(is_bananas)),    median(symm(is_bananas)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'y'); 
plot(median(circ(is_compact)),    median(symm(is_compact)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot(median(circ(is_pancake)),    median(symm(is_pancake)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
plot(median(circ(is_multilobe)),   median(symm(is_multilobe)),'ksq', 'MarkerSize', 12, 'MarkerFaceColor', 'm');


legend('bananas', 'compact', 'pancake','multilobe', 'location', 'EastOutside');
title(['user classification vs shape parameters']);
xlabel('circularity');
ylabel('symmetry');

sel_autobananas=find(circ<lim_circ & symm<lim_symm);
sel_autobananas_usr=find(circ<lim_circ & symm<lim_symm & usr~= 0);
sel_autobananas_usr_true=find(circ<lim_circ & symm<lim_symm & usr== 1);
sel_autobananas_usr_false=find(circ<lim_circ & symm<lim_symm & usr==-1);

auto_results.N_cells=length(circ);
auto_results.N_banana=length(sel_autobananas);
auto_results.N_nonbanana=auto_results.N_cells-auto_results.N_banana;

user_results.N_cells=length(shape_data);
user_results.N_banana=length(circ(is_bananas));
user_results.N_compact=length(circ(is_compact));;
user_results.N_pancake=length(circ(is_pancake));;
user_results.N_multilobe=length(circ(is_multilobe));

save([init.datapath_out, init.exp_label, '__AT030_results_proc.mat'], 'user_results', 'auto_results');

%build excel table
data_banana=[circ(is_bananas) symm(is_bananas)];
data_compact=[circ(is_compact) symm(is_compact)];
data_pancake=[circ(is_pancake) symm(is_pancake)];
data_multilobe=[circ(is_multilobe) symm(is_multilobe)];

header=[{'circularity'},{'symmetry'}];
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], header, 'crescent', 'A1');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], data_banana, 'crescent', 'A2');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], header, 'compact', 'A1');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], data_compact, 'compact', 'A2');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], header, 'dispersed', 'A1');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], data_pancake, 'dispersed', 'A2');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], header, 'multilobed', 'A1');
xlswrite([init.datapath_out, init.exp_label, '_AT030_results_proc.xls'], data_multilobe, 'multilobed', 'A2');


%post_process
disp('auto results:');
disp(['number of cells: ', num2str(auto_results.N_cells)]) ;
disp(['bananas, auto: ' , num2str(auto_results.N_banana),'(' ,num2str(round(100*auto_results.N_banana/auto_results.N_cells)) ,'%)']); 
%1 user categories
disp('user judgment results:');
disp(['number of cells: ', num2str(user_results.N_cells)]) ;
disp(['bananas: ' , num2str(user_results.N_banana),'(' ,num2str(round(100*user_results.N_banana/user_results.N_cells)) ,'%)']); 
disp(['compacts: ', num2str(user_results.N_compact),'(' ,num2str(round(100*user_results.N_compact/user_results.N_cells)) ,'%)']); 
disp(['pancakes: ', num2str(user_results.N_pancake),'(' ,num2str(round(100*user_results.N_pancake/user_results.N_cells)) ,'%)']); 
disp(['multilobes: ',num2str(user_results.N_multilobe),'(' ,num2str(round(100*user_results.N_multilobe/user_results.N_cells)) ,'%)']); 
