function collect_final(users)
        %% 2 run per experiment
        %2023_10_09 JK: revisiting code, adding error intervals using me and Tismas
        %clicks. 
        % tiff-source-conditions:
        % image 1=BSG217_002xy1c2_cmle-1 = -TEV tag, +Xylose: 
        % image 3=BSG219_noXyl_001xy1c2_cmle-1 = +TEV tag, -Xylose
        % image 2=BSG219_45min_Xyl_001xy1c2_cmle-1 = +TEV tag, +Xylose, 45 mins
        % image 4=BSG219_90min_Xyl_001xy1c2_cmle-1 = +TEV tag, +Xylose, 90 mins
        %BlobLegend=([{'1.donut'},{'2.crescent'},{'3.compact'},{'4.other'}]);
        runID=['_run_', datestr(now,'yyyy_mm_dd_HH_MM')];
        col_conditions=[ {'-TEV tag, +Xylose'}; 
                         {'+TEV tag, -Xylose'};
                         {'+TEV tag, +Xylose, 45 mins'};
                         {'+TEV tag, +Xylose, 90 mins'};]

        col_image_index=[1 3 2 4]';
        colz=[{'r'},{'b'}];
        
        %'c=crescent'
        %'t=toroid (donut)'
        %per user, condition:
        %     image index
        %     count_t_1 
        %     count_c_1 
        %     tc_ratio_1 
        Nboot_reps=25;
        N_conditions=length(col_conditions);
        %set up export data
        Wrapup_originalcounts_all_users=zeros(N_conditions,3,length(users));
        Toroid_count_wrapup_all_users=zeros(Nboot_reps,N_conditions,length(users));
        Crescent_count_wrapup_all_users=zeros(Nboot_reps,N_conditions,length(users));
        TC_ratio_wrapup_all_users=zeros(Nboot_reps,N_conditions,length(users));
        TC_ratio_scatax=0*TC_ratio_wrapup_all_users;
        
        figure( 'visible', 'on');
        for uzi=1:length(users) %for every user
               for ci=1:4  %for all conditions (4)
                exp=char(col_conditions(ci));
                image_index=col_image_index(ci);
                all_image_IDs_this_user=[];
                all_blob_IDs_this_user=[];
                usr=char(users{uzi});
                cd(usr);
                click_list=dir('blobstalysis*.mat');
                allblobclass=[];
                %collect all clicked files for this user:
                for ii =1:length(click_list);    
                        load(click_list(ii).name,'typecodes','blobclass','headers');
                        allblobclass=[allblobclass ;...
                                        blobclass];
                        all_image_IDs_this_user=[all_image_IDs_this_user ; allblobclass(:,3);];
                        all_blob_IDs_this_user=[all_blob_IDs_this_user ; allblobclass(:,2)];
                end
                %per user, condition:
                %     image index
                %     count_toroid_1 
                %     count_crescent_1 
                %     tc_ratio_1 
                sel=find(all_image_IDs_this_user==image_index);
                IDs_this_condition=all_blob_IDs_this_user(sel);
                binax = 1:4; 
                [hist_counts, h_Nboot, hist_perc, hPERC_err]=bootstrap_user_classes(IDs_this_condition,binax,Nboot_reps);
                
                %bootstrapped counts of toroids, crescents:
                N_toroids_this_condition_Nboot=h_Nboot(1,:);
                N_crescents_this_condition_Nboot=h_Nboot(2,:);

                %counts and ratios:
                Wrapup_originalcounts_all_users(ci,:,uzi)=[hist_counts(1) hist_counts(2) hist_counts(1)/hist_counts(2)]';
                %bootstrapped counts and ratios:    
                Toroid_count_wrapup_all_users(:,ci,uzi)= N_toroids_this_condition_Nboot'; 
                Crescent_count_wrapup_all_users(:,ci,uzi)= N_crescents_this_condition_Nboot'; 
                Ratio_this_condition_Nboot=N_toroids_this_condition_Nboot./N_crescents_this_condition_Nboot;
                TC_ratio_wrapup_all_users(:,ci,uzi)= Ratio_this_condition_Nboot';                
                TC_ratio_scatax(:,ci,uzi)=ci+0.05*randn(Nboot_reps,1);
                cd ..;
               end
               boxplot(TC_ratio_wrapup_all_users(:,:,uzi),'Labels',col_conditions,'Colors', colz{uzi});
               hold on;
               plot(TC_ratio_scatax(:,:,uzi),TC_ratio_wrapup_all_users(:,:,uzi),'ko', 'MarkerSize', 2);
               
        end
        %% collect fig 4e data:
        ColNames=[{'ratio error'}];
        for ci=1:4  %for all conditions (4)
            sheetname=char(col_conditions(ci));
            SheetData0=0;
            SheetData1=(1:Nboot_reps)';
            for uzi=1:length(users) %for every user
                SheetData0=[SheetData0 (Wrapup_originalcounts_all_users(ci,:,uzi))];
                SheetData1= [SheetData1 Toroid_count_wrapup_all_users(:,ci,uzi) ...
                                        Crescent_count_wrapup_all_users(:,ci,uzi) ...
                                        TC_ratio_wrapup_all_users(:,ci,uzi)];
            end
            ColNames=[{'resample run'} , {'toroid counts user 1'},  {'crescent counts user 1'}, {'ratio user 1'},...
                                         {'toroid counts user 2'},  {'crescent counts user 2'}, {'ratio user 2'}];
                                     
            xlswrite(['blobstalyzer_all_users_wrapup_final',runID,'.xlsx'], ColNames,sheetname,'A1');
            xlswrite(['blobstalyzer_all_users_wrapup_final',runID,'.xlsx'], SheetData0,sheetname,'A2');   
            xlswrite(['blobstalyzer_all_users_wrapup_final',runID,'.xlsx'], [{'original'}],sheetname,'A2');
            xlswrite(['blobstalyzer_all_users_wrapup_final',runID,'.xlsx'], SheetData1,sheetname,'A3');   
        end                     
         saveas(gcf,['blobstalyzer_all_users__wrapup',runID,'.jpg']);
         plot2svg(['blobstalyzer_all_users__wrapup',runID,'.svg'], gcf);