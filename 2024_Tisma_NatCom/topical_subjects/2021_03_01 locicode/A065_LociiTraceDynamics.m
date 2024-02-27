function A065_LociiTraceDynamics(initval)
%JWJK_A:-------------------------------------------------------------------
%
%Building of  locii traces
%
%Summary: code loads image stack per cropped cell and performs spot
%analysis on each.
%
%Approach: first, stack is loaded of rfp or cfp. average of stack is used
%for defining a fixed, smaller ROI. This is then used for 2D-Xcor tracking
%
%Input: movie stack onf cropped cell images; all channels 
%
%Output: collection of traces per cell
%
%References: written by JWJK
%:JWJK_A-------------------------------------------------------------------

close all
cur_path=pwd;
pix2nm=initval.pix2nm;
Labels=initval.SaveLabels(2:4);
SimbolColors=[{'r'},{'b'}, {'k'}];
%Labels=[{'rfp'}];
LB=length(Labels);
txtindir=strcat(initval.imagepath,'A60_LociiTraces\');
for Lbi=1:LB                        %for all label types
    ThisLabel=char(Labels(Lbi)); 
    cd(txtindir);
    tracelist=dir(['Trace*', ThisLabel,'.txt']);
    cd(cur_path);
    LTS=length(tracelist);   
    for ii=1:LTS %for all traces
        disp(strcat(initval.expname,':',initval.subdir,'-','_Label_',ThisLabel,'_MSDofTrace',num2str(ii),'of',num2str(LTS)));
        xytrace=pix2nm*dlmread(strcat(txtindir,'\',tracelist(ii).name));
        if ii==1;
            [LT,~]=size(xytrace); 
            MSDs=NaN*zeros(LTS,5*LT);   %no of traces _ lag time 
            %note: (create 5xmargin for long traces); NaN will not
            %contribute
            AllTracesX=NaN*zeros(5*LT,LTS);
            AllTracesY=NaN*zeros(5*LT,LTS);
        end     
       AllTracesX(1:LT,ii)=xytrace(:,1); %each column is a trace
       AllTracesY(1:LT,ii)=xytrace(:,2);  
       
       [LT,~]=size(xytrace); 
       for jj=1:LT  %for all available time differences in this trace
            MSDs(ii,jj)=Get_MeanSquareDisplacement(xytrace(:,1),xytrace(:,2),jj);
            %each row is MSD output of one race
       end    
    end
       
    MSD_Av=nanmean(MSDs);
    MSD_Std=nanstd(MSDs);
    
    %crop empty upper parts
    sel=find(~isnan(MSD_Av));     
    AllTracesX=AllTracesX(sel,:);
    AllTracesY=AllTracesY(sel,:);
    MSDs=MSDs(:,sel);
    MSD_Av=(MSD_Av(sel));
    MSD_Std=(MSD_Std(sel));    
    axz=(1:length(sel));


    %% plot menu
   clr=char(SimbolColors(Lbi));
    figure(Lbi);
    subplot(2,2,1);
        plot(axz,AllTracesX,strcat(clr,'-')), hold on;
        title(ThisLabel);
        ylabel('position, nm')
        xlabel('frame index');
        legend('X-traces');    

    subplot(2,2,2);
        plot(axz,AllTracesY,strcat(clr,'-')), hold on;
        title('time traces')
        ylabel('position, nm')
        xlabel('frame index');
        legend('Y-traces');

    subplot(2,2,3);
        plot(AllTracesX,AllTracesY,strcat(clr,'o'), 'MarkerSize',2), hold on;
        axis equal;
        title('all positions');   
        ylabel('position, nm')
        xlabel('position, nm');
        pause(0.01);

    subplot(2,2,4);    
        plot(axz',MSDs','k-'), hold on;
        plot(axz',MSD_Av',strcat(clr,'o-'), 'MarkerSize',7), hold on;
        title('Means square displacement');   
        xlabel('lag time, frame units')
        ylabel('msd, nm^2');
        %ylim([0 0.5E5]);

    outdata=[axz' MSD_Av' MSD_Std' MSDs'];    
    %% saving
    ColNames=[{'lag time'} , {'average'}, {'standard deviation'}, repmat({'single locus curve'},1,LTS)];
    xlswrite(strcat(initval.imagepath,'Exp_dynamics','_A065_MSD_Curves_',ThisLabel,'.xlsx'),ColNames,'Sheet1','A1');
    xlswrite(strcat(initval.imagepath,'Exp_dynamics','_A065_MSD_Curves_',ThisLabel,'.xlsx'),outdata,'Sheet1','A2');
    saveas(gcf, strcat(initval.imagepath,'Exp_dynamics','_A065_MSD_Curves_',ThisLabel,'.jpg'));
end

