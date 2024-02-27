function collect_rebuttal(users)
        %% 2 run per experiment
        %2023_10_09 JK: revisiting code, adding error intervals using me and Tismas
        %clicks. 
        % tiff-source-conditions:
        % image 1=BSG217_002xy1c2_cmle-1 = -TEV tag, +Xylose: 
        % image 3=BSG219_noXyl_001xy1c2_cmle-1 = +TEV tag, -Xylose
        % image 2=BSG219_45min_Xyl_001xy1c2_cmle-1 = +TEV tag, +Xylose, 45 mins
        % image 4=BSG219_90min_Xyl_001xy1c2_cmle-1 = +TEV tag, +Xylose, 90 mins
        %BlobLegend=([{'1.donut'},{'2.crescent'},{'3.compact'},{'4.other'}]);
        col_conditions=[ {'-TEV tag, +Xylose'}; 
                         {'+TEV tag, -Xylose'};
                         {'+TEV tag, +Xylose, 45 mins'};
                         {'+TEV tag, +Xylose, 90 mins'};]

        col_image_index=[1 3 2 4]';

        
        %'c=crescent'
        %'t=toroid (donut)'
        %per user, condition:
        %     image index
        %     count_t_1 
        %     count_c_1 
        %     tc_ratio_1 
        runID=['_run_', datestr(now,'yyyy_mm_dd_HH_MM')];
        wrapup_all_users=zeros(4,4,length(users));
        for uzi=1:length(users) %for every user
               wrapup_this_user=zeros(4,4);
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
                N_toroids_this_condition=length(find(IDs_this_condition==1));
                N_crescents_this_condition=length(find(IDs_this_condition==2));
                Ratio=N_toroids_this_condition/N_crescents_this_condition;
                wrapup_this_user(ci,:)=[image_index ...
                    N_toroids_this_condition...
                    N_crescents_this_condition...
                    Ratio];
                cd ..;
                dum=1;
               end
             wrapup_all_users(:,:,uzi)= wrapup_this_user;
        end
        %% now, compare
        for ci=1:4  %for all conditions (4)
            nme=col_conditions{ci};
            toroid_counts=squeeze(wrapup_all_users(ci,2,:));
            crescent_counts=squeeze(wrapup_all_users(ci,3,:));
            ratios=squeeze(wrapup_all_users(ci,4,:));
            avratios(ci)=mean(ratios);
            err_ratios(ci)=(max(ratios)-min(ratios))/2;
            dum=1;
                                 
        end
        close all;
        figure( 'visible', 'on');
        bar((1:4),avratios); hold on;
        errorbar((1:4),avratios, err_ratios, 'rsq');
        
        wrapup=[col_image_index avratios' err_ratios' ];
             
        ColNames=[{'condition'} , {'tiff_index'}, {'average ratio'} ,{'ratio error'}];
        xlswrite(['blobstalyzer_all_users_wrapup_rebuttal',runID,'.xlsx'], ColNames,'Sheet1','A1');
        xlswrite(['blobstalyzer_all_users_wrapup_rebuttal',runID,'.xlsx'], col_conditions,'Sheet1','A2');
        xlswrite(['blobstalyzer_all_users_wrapup_rebuttal',runID,'.xlsx'],wrapup,'Sheet1','B2');
        
         saveas(gcf,['blobstalyzer_all_users_wrapup_rebuttal',runID,'.jpg']);
         plot2svg(['blobstalyzer_all_users_wrapup_rebuttal',runID,'.svg'], gcf);