function A400_Correlation
%JWJK_A:-------------------------------------------------------------------
%Title: spatial and temporal correlation analysis
%
%Summary: % Spatial autocorrelation analysis was performed on 10 individual images per movie. 
% For each autocorrelation output image, a radial average was recorded 
% starting from the main central correlation peak. 
% The resulting spatial radial correlation curve was subjected to 
% maxima analysis. The first maximum after radius R=0  indicated 
% the most predominant distance between wave edges, 
% irrespective of propagation direction. 
% This distance was denoted as lambda. For temporal correlation, 
% we generated 20 x-t or y-t kymographs per movie 
% (10 in ‘x’ direction and 10 in ‘y’ direction) 
% evenly distributed over the middle 0.7 fraction of an image. 
% For each such kymograph, an autocorrelation analysis was performed. 
% The x=0 or y=0 line of these x-t autocorrelation maps then in effect 
% represents a temporal correlation curve averaged over all the original 
% image points on this line. Next, these correlation curves were 
% median averaged between different kymographs. 
% Thus, the final correlation curve in effect represents the average 
% temporal correlation signal sampled from 20x512 surface locations. 
% Analogous to the spatial correlation analysis, the first maximum after t=0 
% indicated a main oscillation period. 
% Experimental repeats of the same conditions (concentrations and height) 
% were median averaged.            
%
%Input: directory with pre-cleaned movies and directory with tilted xt or yt 
%versions of these movies, named following experimental
%conditions. as follows:
% First letter of the filename means height, ABCDEF from lowest to highest.
% translates as heights: [2 6 8 15 25 57] microns
% second comes number which will be E/D ratio 1,2,3,4,5 , 
% translates as ratios: (0.5, 0.75. 1, 2, 3) 
% Last letter is channel (D -MinD, E -MinE)
% Second last number is experiment session
% Last number is repeat of experiment session with same conditions.
% 
% Example: 'C_3_E_1_2' means: 
%     height 8, ratio 1,MinE signal, session 1, repeat 2.

%Output: 
% 1) C_3_E_1_2_cln_overview_Height08Ratio2.00ChanERepeat01.jpg: plots
% summarizing the results of the spatial and temporal periodicity analysis

% 2) A400_Correlation_results_permovie.mat:
% Data output description of Greg's 'correlation' data
%  % General -----------------------------------------
%  MinED_mov = 1×312 struct array with fields:
%      filname: name of movie
%      channel: colorlabel  channel ('D' or 'E')
%      height: sampe height
%      EDratio: MinE-MinD ratio
%      repeat: experimental repeat index of this EDratio-height combination
%      xy: see below
%      xt: see below
%      yt: see below
%      info: 
%          frames2mins: 1
%           pix2nm: 594
%   
% Details ----------------------------------------------
% For example, MinED_mov(1).xy  contains data of the first movie.      
%  Each field xy or xt contains, per movie, some details of the obained correlation curves.
% The xy correlation curves are obtained from radial averages of the
% correlation maps, the xt correlation curves are taken along the x=0 line
% (since then, we look at the autocorrelation of a wave with itself)
%          suffix _mv means per movie frame (here, 10 were chosen)
%          suffix _av means average of these values
% data:
%          decay_val_labda_av:
%          decaylength_abs_av: 
%          decaylength_sign_av: 
%          decayperiodicityfactor_av: 
%          envelope_abs: empty
%          example: [1×1 struct]: 
%          mn1pos_av: position of first minimum (''anticorrelation'')
%          mn1val_av: amplitude of this peak
%          pk2pos_av: position of first maximum (indicative of periodicity)
%          pk2val_av: amplitude of this peak
%          r_axis: for plotting: axis values
%          radialprofiles: for plotting: profiles
%          strengthofmodulation_av: factor relating the modulation depth of the
%          correlation peaks to the overall decay value, to express how
%          'periodic' the waves look
 
% 3) A400_Correlation_results_tabulated.mat: here, the above parameters are
% saved as columns, one row per movie

% 4) A400_Correlation_results_tabulated.xlsx: as the .mat file

% 5) A400_Correlation_results_matrices..jpg: visual representation of
% parameters, set in a matrix of heights and MinE/MinD ratio.


%Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
%code designed & written Jacob Kerssemakers 2016 
%:JWJK_A-------------------------------------------------------------------




addpath(genpath(pwd));
close all;



actions.test='test';
actions.test='development';

actions.movie_process_it =1;
actions.tabulate_it=1;
actions.plot_it=0;  
workpath=pwd;

matrixoptions.flattenregimes=0;


initval.pix2nm=594;

initval.kymoband=0.8; %pick 8 kymographs from middle 0.7 part of image
switch actions.test
    case 'Sabrina_1'
        initval.maindirname=[workpath,'\testdata\'];
        initval.experimentseriesname='movies_test';  
        initval.dorepeats=1;
        initval.reps_permovie=1;
        initval.reps_per_kymostack=10;
        set(figure(1), 'visible','on');
        heightlabels=[{'C'}];
        heights=[8];
        ED_ratiolabels=[{'3'}];
        EDratios=[2];
        EDfr2mins=[1];
        Channels=[{'E'}];
end


outdirname=strcat(initval.maindirname,initval.experimentseriesname,'_matlaboutput\');
file_index=0;  %file counter

if actions.movie_process_it 
    tic
    if ~isdir(outdirname), mkdir(outdirname); end
    %output columns to export
    MinED_mov=struct(...
    'filname' ,[],...
    'channel' ,[],... 
    'height', [],...
    'EDratio',[],...
    'repeat',[],...
    'xy',[],... %properties related to xy-images
    'xt',[],... %properties related to xt-kymographs
    'yt',[]...  %properties related to yt-kymographs
    );


    %% work all movies
    for ch=1:length(Channels)               %for all channels
        CH=char(Channels(ch));
        for hg=1:length(heightlabels)       %for all heights        
            HG=char(heightlabels(hg));
            hght=heights(hg);
            for rt=1:length(ED_ratiolabels) %for all ratios
                RT=char(ED_ratiolabels(rt)); %(for now, just one image)
                EDratio=EDratios(rt);
                fr2min=EDfr2mins(rt);
                %find and load             
                filrepname=strcat(HG,'_',RT,'_',CH,'*cln.tif');
                %this is the search template
                reslicefilname=strcat(HG,'_',RT,'_',CH,'*cln_resliced.tif');           
                
                cd(strcat(initval.maindirname,initval.experimentseriesname, '_cln\'));
                filelist=dir(filrepname); %list; including repeats of same measurement!
                cd(pwd);
                if ~isempty(filelist); 
                if initval.dorepeats, 
                    Reps=length(filelist); 
                else
                    Reps=1;
                end
                
                    for rp=1:Reps
                    disp(strcat(...
                    num2str(length(Channels)-ch+1), 'channels',...
                    num2str(length(heightlabels)-hg+1),'heights',...
                    num2str(length(ED_ratiolabels)-rt+1),'ratios',...
                    num2str(Reps-rp+1), 'repeats',...
                            'to go')); 
                    filname=filelist(rp).name                   
                    file_index=file_index+1;  %counter
                    MinED_mov(file_index).filname=filname;
                    MinED_mov(file_index).channel=CH; 
                    MinED_mov(file_index).height=hght;
                    MinED_mov(file_index).EDratio=EDratio;
                    MinED_mov(file_index).repeat=rp;
                    MinED_mov(file_index).info.frames2mins=fr2min;
                    MinED_mov(file_index).info.pix2nm=initval.pix2nm;

                    %% analyze I: regular movie; first  ... images
                    ThisMov_XY=Get_XY_Movie_properties(filname,initval);
                    MinED_mov(file_index).xy=ThisMov_XY;

                    %% analyze II: resliced (tilted) movies (made before)
                    ThisKymo_XT=Analyse_Kymo_Band(initval,filname,'x',fr2min,rp);
                    MinED_mov(file_index).xt=ThisKymo_XT;              
                    ThisKymo_YT=Analyse_Kymo_Band(initval,filname,'y',fr2min,rp);
                    MinED_mov(file_index).yt=ThisKymo_YT;
                    end
                end
            end
        end
    end
    toc
    save(strcat(outdirname,'A400_Correlation_results_permovie.mat'),'MinED_mov');  
end              
                
                
 if actions.tabulate_it|actions.plot_it  
    load(strcat(outdirname,'A400_Correlation_results_permovie.mat'),'MinED_mov'); 
    [~,fils]=size(MinED_mov);
 end
 
 cd(workpath);
    
 if actions.tabulate_it
    for fl=1:fils  %build columns for excel table (or plotting use)
        ThisMov=MinED_mov(fl);
        filnameout(fl)={ThisMov.filname};
        chanout(fl)={ThisMov.channel};
        hgt_out(fl)=ThisMov.height;
        EDratio_out(fl)=ThisMov.EDratio;
    
        %spatial correlation
        xycor_PK1radii(fl)=ThisMov.xy.cor.pk2pos_av;
        xycor_PK1vals(fl)=ThisMov.xy.cor.pk2val_av;
        xycor_MN1radii(fl)=ThisMov.xy.cor.mn1pos_av;
        xycor_MN1vals(fl)=ThisMov.xy.cor.mn1val_av;
        
        xycor_decaylength_abs(fl)=ThisMov.xy.cor.decaylength_abs_av;
        xycor_decaylength_sign(fl)=ThisMov.xy.cor.decaylength_sign_av;
        xycor_decayperiodicityfactor(fl)=ThisMov.xy.cor.decayperiodicityfactor_av;
        xycor_strengthofmodulation(fl)=ThisMov.xy.cor.strengthofmodulation_av;
        
        %temporal correlation in x and y direction averaged
        xtcor_PK1radii(fl)=(ThisMov.xt.cor.pk2pos_av+...
                            ThisMov.yt.cor.pk2pos_av)/2;
        xtcor_PK1vals(fl)=(ThisMov.xt.cor.pk2val_av+...
                          ThisMov.yt.cor.pk2val_av)/2;
        xtcor_MN1radii(fl)=(ThisMov.xt.cor.mn1pos_av+...
                            ThisMov.yt.cor.mn1pos_av)/2;
        xtcor_MN1vals(fl)=(ThisMov.xt.cor.mn1val_av+...
                          ThisMov.yt.cor.mn1val_av)/2;
        xtcor_decaylength_abs(fl)=(ThisMov.xt.cor.decaylength_abs_av+...
                                   ThisMov.yt.cor.decaylength_abs_av)/2;
        xtcor_decaylength_sign(fl)=(ThisMov.xt.cor.decaylength_sign_av+...
                                    ThisMov.yt.cor.decaylength_sign_av)/2;
        xtcor_decayperiodicityfactor(fl)=(ThisMov.xt.cor.decayperiodicityfactor_av+...
                                          ThisMov.yt.cor.decayperiodicityfactor_av)/2;
        xtcor_strengthofmodulation(fl)=(ThisMov.xt.cor.strengthofmodulation_av+...
                                        ThisMov.yt.cor.strengthofmodulation_av)/2;
               
              
        %boundary counts
        xtcor_domainL(fl)=(ThisMov.xt.fourier.xtdomains_av+...
                             ThisMov.xt.fourier.xtdomains_av)/2*initval.pix2nm/1000;
        
    end
    %% Saving
    
    %% .MAT  Column Save
    save(strcat(outdirname,'A400_Correlation_results_tabulated.mat'),'heights','EDratios');
    if 0
     save(strcat(outdirname,'A400_Correlation_results_tabulated.mat'),...
                                   'filnameout','chanout','hgt_out', 'EDratio_out',...
                                  'xycor_MN1radii', 'xycor_MN1vals', 'xycor_PK1radii', 'xycor_PK1vals',...
                                  'xycor_decaylength_abs', 'xycor_decaylength_sign',...
                                  'xycor_decayperiodicityfactor', 'xycor_strengthofmodulation',...
                                  'xtcor_MN1radii', 'xtcor_MN1vals', 'xtcor_PK1radii', 'xtcor_PK1vals',...
                                  'xtcor_decaylength_abs', 'xtcor_decaylength_sign',...
                                  'xtcor_decayperiodicityfactor', 'xtcor_strengthofmodulation',...                                 
                                  'xtcor_domainL', '-append');
  
    end
   %% quick matrix style (median of groups)
    
    [xycor_PK1radii_matrix,~]=QuickMatrix_it(xycor_PK1radii,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions);    
    [xtcor_PK1radii_matrix,~]=QuickMatrix_it(xtcor_PK1radii,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions);        
    [xycor_decaylength_abs_matrix,~]=QuickMatrix_it(xycor_decaylength_abs,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions); 
    [xtcor_decaylength_abs_matrix,~]=QuickMatrix_it(xtcor_decaylength_abs,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions);    
    [xycor_decayperiodicityfactor_matrix,~]=QuickMatrix_it(xycor_decayperiodicityfactor,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions); 
    [xtcor_decayperiodicityfactor_matrix,~]=QuickMatrix_it(xtcor_decayperiodicityfactor,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions);    
    [xtcor_domainL_matrix,mask_matrix]=QuickMatrix_it(xtcor_domainL,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions); 
     
    save(strcat(outdirname,'A400_Correlation_results_tabulated.mat'),...
        'xycor_PK1radii_matrix',...
        'xycor_decaylength_abs_matrix',...
        'xycor_decayperiodicityfactor_matrix',...       
        'xtcor_PK1radii_matrix',...
        'xtcor_decaylength_abs_matrix',...
        'xtcor_decayperiodicityfactor_matrix',...
        'xtcor_domainL_matrix',...
        'mask_matrix',...
        '-append');
         
    figure(215);
    subplot(2,4,2); pcolor(xycor_PK1radii_matrix); colormap hot; shading flat; axis off; title('X-labda');
    subplot(2,4,3); pcolor(xycor_decaylength_abs_matrix); colormap hot; shading flat; axis off;title('X-abs.decay');
    subplot(2,4,4); pcolor(xycor_decayperiodicityfactor_matrix); colormap hot; shading flat;axis off; title('X-strength');
    
    subplot(2,4,8); pcolor(xtcor_domainL_matrix); colormap hot; shading flat; axis off;title('Phasedomains');    
    
    subplot(2,4,5); pcolor(xtcor_PK1radii_matrix); colormap hot; shading flat; axis off;title('T-tau');
    subplot(2,4,6); pcolor(xtcor_decaylength_abs_matrix); colormap hot; shading flat;axis off; axis off;title('T-abs.decay');
    subplot(2,4,7); pcolor(xtcor_decayperiodicityfactor_matrix); colormap hot; shading flat; axis off;title('T-strength');
    
    subplot(2,4,1); pcolor(mask_matrix); colormap hot; shading flat; axis off; title('Regimes');
 
    %save
    figoutname=strcat(outdirname,'A400_Correlation_results_matrices.jpg');
    saveas(gcf,figoutname, 'jpg');
    pause(0.1);
    
    
                              
    %% EXCELSAVE
    %save main results in columns for easy excel handling
    ExcelColNames=...
    [{'filename'}, {'channel'}, {'height,mu'} ,{'ED-ratio'},...
    {'xycor mn1pos,mu'},{'xycor mn1val'},...
    {'xycor peak1pos,mu'},{'xycor peak1val'},{'xycor peakvalley-diff'},...
    {'xtcor mn1pos,mins'},{'xtcor mn1val'},...
    {'xtcor peak1pos,mins'},{'xtcor peak1val'},{'xtcor peakvalley-diff'},...
    {'ytcor mn1pos,mins'},{'ytcor mn1val'},...
    {'ytcor peak1pos,mins'},{'ytcor peak1val'},{'ytcor peakvalley-diff'}...
    ];
    
          
    xlswrite(strcat(outdirname,'A400_Correlation_results_tabulated.xlsx'),ExcelColNames,'Sheet1','A1');
    xlswrite(strcat(outdirname,'A400_Correlation_results_tabulated.xlsx'),filnameout','Sheet1','A2');
    xlswrite(strcat(outdirname,'A400_Correlation_results_tabulated.xlsx'),chanout','Sheet1','B2');
    xlswrite(strcat(outdirname,'A400_Correlation_results_tabulated.xlsx'),...
        [hgt_out' EDratio_out' ...
        xycor_MN1radii' xycor_MN1vals' xycor_PK1radii' xycor_PK1vals' xycor_PK1vals'-xycor_MN1vals'...
        xtcor_MN1radii' xtcor_MN1vals' xtcor_PK1radii' xtcor_PK1vals' xtcor_PK1vals'-xtcor_MN1vals'],...
        'Sheet1','C2');    
 end
     
 if actions.plot_it
     close all;
        for fl=1:fils  %build columns for excel table (or plotting use)
        
        ThisMov=MinED_mov(fl);
        filename=ThisMov.filname;
        chan=ThisMov.channel;
        hght=ThisMov.height;
        EDratio=ThisMov.EDratio;
        rep=ThisMov.repeat;
        

        %show X-X
        titl=strcat('Height',num2str(hght,'%02.0f'),'Ratio',num2str(EDratio,'%01.2f'),'Chan',chan,'Repeat',num2str(rep,'%02.0f'));
        subplot(2,3,1); pcolor(ThisMov.xy.cor.example.plotim); shading flat; colormap hot; 
        axis equal; axis tight;
        title(titl); xlabel('position'); ylabel('position');
        subplot(2,3,2); pcolor(ThisMov.xy.cor.example.cormap); shading flat; colormap hot;
        axis equal; axis tight;
        title('autocor X-X');

        subplot(2,3,3); 
%          plot(ThisMov.xy.cor.r_axis',...
%               ThisMov.xy.cor.envelope_abs','k-','Linewidth',1); hold on;
        plot(ThisMov.xy.cor.r_axis',...
             ThisMov.xy.cor.radialprofiles','Linewidth',1); hold off;
        
        title('autocor');
%         text(40,0.5,strcat('R:',num2str(ThisMov.xy.cor.pk2rad_av)));
%         text(40,0.3,strcat('C_v_a_l:',num2str(ThisMov.xy.cor.pk2val_av)));
        axis tight
        ylim([-0.5 1]);
        xlabel('distance, pixels');
        ylabel('C, a.u.');


        %X-t
        subplot(2,3,4); pcolor(ThisMov.xt.cor.example.plotim); shading flat; colormap hot;
        xlabel('frame'); ylabel('position');
        axis square, axis tight;            
        title('X-T plot');

        subplot(2,3,5); pcolor(ThisMov.xt.cor.example.cormap_xt); shading flat; colormap hot;
        axis square, axis tight;
        %title('autocor X-T');

        subplot(2,3,6); 
%         plot(ThisMov.xt.cor.t_axis',...
%              ThisMov.xt.cor.envelope_abs','k-','Linewidth',1); hold on;
        plot(ThisMov.xt.cor.t_axis',...
             ThisMov.xt.cor.prf_tcor','Linewidth',1); hold off;

%          text(10,0.4,strcat('R:',num2str(ThisMov.xt.cor.pk2pos_av)));
%          text(10,0.2,strcat('C_v_a_l:',num2str(ThisMov.xt.cor.pk2val_av)));
        axis tight;
        xlabel('time, minutes');
        ylabel('C, a.u.');

        %save
        matrixname=titl;
        figoutname=strcat(outdirname,filename(1:end-4),'_overview_',matrixname,'.jpg');
        saveas(gcf,figoutname, 'jpg');
        pause(0.1);
        end


end






function stetsoncurve=BuildStetsonCurve(radialprofile);
%Build a flat-topped edge mask
LL=length(radialprofile);
dropofflengthIn=ceil(LL*0.2);
dropofflengthOut=ceil(LL*0.2);
flatlength=LL-dropofflengthIn-dropofflengthOut;
EdgesIn=(hann(2*dropofflengthIn))';
EdgesOut=(hann(2*dropofflengthOut))';
stetsoncurve=[...
    EdgesIn(1:dropofflengthIn) ...
    ones(1,flatlength)...
    EdgesOut(dropofflengthOut+1:end)];
dum=1;

function [crmx, r_axis,radialprofile,radialprops]=RadialAutoCor(im);
%find spatial autocorrelation
maxsearchrange=[10  200];
 
%1 normalize
    im=im-min(im(:));
    im=im/sum(im(:));
    im=im-mean(im(:));

    mirror_im=fliplr(flipud(im));    
    crmx=fftshift(real(ifft2(fft2(im).*conj(fft2(im)))));      %auto-correlation
    crmx=crmx/max(crmx(:));
 %sample from center
    QI=TrackXY_by_QI_Init(crmx); 

         [rr,cc]=size(crmx);
         hfzy=floor(rr/2);
         hfzx=floor(cc/2);
         x0=hfzx+1;
         y0=hfzy+1;
         Xsamplinggrid=QI.X0samplinggrid+x0;
         Ysamplinggrid=QI.Y0samplinggrid+y0;
         allprofiles=(interp2(crmx,Xsamplinggrid,Ysamplinggrid,'NaN'));

         %translate axes
         r_axis=QI.radbii;
         radialprofile_sum=nansum(allprofiles);
         countprofile=sum(~isnan(allprofiles)*1.0);
         radialprofile=radialprofile_sum./countprofile;  %weigths
         radialprofile=radialprofile/nanmax(radialprofile);  
%          if 1
%              figure(2); 
%              plot(radialprofile);
%              dum=1;[~]=ginput(1);
%          end
         %check for maximum
         smz=10;
         profileprops=AnalyzeCorProfile(r_axis,radialprofile,smz);         
         radialprops=profileprops;
        
      
 function [crmx,t_axis,prf_tcor,profileprops]=XT_AutoCor(xt_im);     
%find spatial autocorrelation
%to add: ''central cross'' repeats


%to ad: once main period is found, get phase of this line
%get a histogram of this phase
%if the histogram is 'flat', all phases are present equally, which is
%representative of a running wave
%if it's [very] discrete, the area (or line) behaves in areas of 'single pahse''
 
%1 normalize
    xt_im=xt_im-min(xt_im(:));
    xt_im=xt_im/sum(xt_im(:));
    xt_im=xt_im-mean(xt_im(:));

    mirror_im=fliplr(flipud(xt_im));    
    crmx=(real(ifft2(fft2(xt_im).*conj(fft2(xt_im)))));      %auto-correlation
    [rr,cc]=size(crmx);   
    prf_tcor=squeeze(crmx(1,1:round(cc/2)));
    t_axis=1:length(prf_tcor);
    prf_tcor=prf_tcor/max(prf_tcor);
    smz=2;
    profileprops=AnalyzeCorProfile(t_axis,prf_tcor,smz);
    crmx=fftshift(crmx);

    
    
  function profileprops=AnalyzeCorProfile(r_axis,radialprofile_ori,smz)
         %step 0: clean up for fit
         sel=find(~isnan(radialprofile_ori));
         radialprofile_ori=radialprofile_ori(sel);
         r_axis=r_axis(sel);
         sel1=find(~isnan(r_axis));
         radialprofile_ori=radialprofile_ori(sel1);
         r_axis=r_axis(sel1);
         
         if ~isempty(radialprofile_ori)
         
         %step 1: measure absolute average decay         
         radialprofile_0shift=radialprofile_ori-radialprofile_ori(end);
         radialprofile_abs=abs(radialprofile_0shift);
         
         radialprofile_abs=radialprofile_abs/radialprofile_abs(1);
         radialprofile_0shift=radialprofile_0shift/radialprofile_0shift(1);
         
         decayfittype = fittype('(exp(-x/delta))', ...
             'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'delta'});
         %fit on absolute correlation
         startpnt=sum(r_axis.*radialprofile_abs)/sum(radialprofile_abs);
         prms_abs = fit(r_axis(1:end-1)',radialprofile_abs(1:end-1)',decayfittype,...
             'Lower',[1],'Upper',[max(r_axis')],'StartPoint',startpnt,...
             'MaxFunEvals',250);
         decayfit_abs=exp(-r_axis/prms_abs.delta);
         %fit on relative correlation
         prms_signed = fit(r_axis(1:end-1)',radialprofile_0shift(1:end-1)',decayfittype,...
             'Lower',[1],'Upper',[max(r_axis')],'StartPoint',startpnt,...
             'MaxFunEvals',250);
         decayfit_signed=exp(-r_axis/prms_signed.delta);
         
         
         %step 2: assuming periodicity, measure position of the first
         %valleys and peaks
         radialprofile=radialprofile_ori;
         %check for maximum   
         %1 find minima and maxima
         goodpeaks=1;        
         rpm1=smooth(radialprofile(1:end-2),smz);
         rp0=smooth(radialprofile(2:end-1),smz);
         rpp1=smooth(radialprofile(3:end),smz);
         ax0=r_axis(2:end-1);
         
         selmaxima=find((rpm1<rp0)&(rpp1<rp0));
         axmax=ax0(selmaxima); %associated radial positions                
         selminima=find((rpm1>rp0)&(rpp1>rp0));
         if ~isempty(selminima)
         idx2=find(selmaxima>selminima(1));
         
         mn1val=rp0(selminima(1));
         mn1pos=ax0(selminima(1));
         
             if~isempty(idx2)  %are there maxima beyond first minimum?
             selmaxima_2=selmaxima(idx2);
             pks=rp0(selmaxima_2);
             pks_ax=ax0(selmaxima_2);


             [pk2val,idx3]=max(pks); 
             pk2pos=pks_ax(idx3);  
             %maximum peak beyond 1st minimum

             sel=find(ax0==pk2pos);   
             if (sel>1)&(sel<length(ax0)&~isempty(sel)) %subpixel step
                pksection=rp0(sel-1:sel+1);       
                x=subpix_aroundzero(pksection); %correction    
                pksax=ax0(sel-1:sel+1);
                pk2pos=interp1(-1:1,pksax,x);
             end
             else  %off limits, put at twice anti-correlation
                 goodpeaks=1;                
                 pk2pos=2*mn1pos;
                 pk2val=NaN;
             end 
         else
             goodpeaks=0;
         end
         if goodpeaks==0
             pk2val=NaN; pk2pos=NaN;  %no secondary maximum
             mn1val=NaN; mn1pos=NaN;
         end
         
         profileprops.decaylength_abs=prms_abs.delta;
         profileprops.decaylength_sign=prms_signed.delta;
         
         profileprops.mn1pos=mn1pos;
         profileprops.mn1val=mn1val;
         profileprops.pk2pos=pk2pos;
         profileprops.pk2val=pk2val;
         
         %Now, define a modulation depth
         profileprops.strengthofmodulation=profileprops.pk2val-profileprops.mn1val;
         %Compare this with the value of the absolute decay fit between
         %these two extrema:
         profileprops.decay_val_labda=(exp(-pk2pos/prms_abs.delta)+exp(-mn1pos/prms_abs.delta))/2;
         %if the profile strongly modulates and has no non-periodic
         %component, this decay value is approximately half of the
         %modulation depth
         profileprops.decayperiodicityfactor=profileprops.strengthofmodulation/profileprops.decay_val_labda;
         
         profileprops=orderfields(profileprops);
         
         else
         profileprops.decaylength_abs=NaN;
         profileprops.decaylength_sign=NaN;
         profileprops.decayperiodicityfactor=NaN;
         profileprops.strengthofmodulation=NaN;
         profileprops.decay_val_labda=NaN;
         
         profileprops.mn1pos=NaN;
         profileprops.mn1val=NaN;
         profileprops.pk2pos=NaN;
         profileprops.pk2val=NaN;
             
         end

         if 0
             close all
%              subplot(1,2,1);
%              plot(r_axis,radialprofile_0shift,'b-'); hold on;
%              plot(r_axis,decayfit_signed,'k-'); hold off;
%              xlabel('position, a.u.');
%              ylabel('correlation value, a.u.');
%              legend('correlation', 'exp. fit');
%              axis tight
%              subplot(1,2,2);
             plot(r_axis,radialprofile_0shift,'b-'); hold on
             plot(r_axis,radialprofile_abs,'k-'); hold on
             plot(r_axis,decayfit_abs,'r-'); hold on
             legend('correlation','absolute correlation', 'exp. fit');
             xlabel('position, a.u.');
             ylabel('correlation value, a.u.');
             axis tight
             profileprops
             [~]=ginput(1);
             
             dum=1;
         end
          
    
    function  x=subpix_aroundzero(prfx);
         %3-point subpixel fit
         xax=[-1:1:1]'; [~,mxi]=max(prfx);
         lpr=length(prfx);
         idxes=mxi-1:1:mxi+1;
         sel=find(idxes<1); idxes(sel)=idxes(sel)+lpr;
         sel=find(idxes>lpr); idxes(sel)=idxes(sel)-lpr;
         prfx=prfx(idxes);   %peak parabols with edge transfer     
         prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
    
    function ThisKymo_XT=Analyse_Kymo_Band(initval,filname,whichaxis,fr2min,rp);
        %We analyze groups of kymographs belonging to this movie, 
        %taken as a band around the middle xy area.
        ThisKymo_XT=struct('kymopositions',[]);                
        %get x and y stack info
        reslicefilname=strcat(filname(1:end-4),'_resliced_',whichaxis,'.tif');    
        slicepth=strcat(initval.maindirname,initval.experimentseriesname,'_cln_rs\',reslicefilname);
        firstslice=double(imread(slicepth,'Index',1));
        if strcmp(whichaxis, 'y'), firstslice=firstslice'; end
            
        info = imfinfo(slicepth);    
        [Ly,~]=size(info);                                           
        [Lt,Lx]=size(firstslice);               
        %define a band of x-cross-sections trhough the image middle
        Nsect=initval.reps_per_kymostack; band=initval.kymoband; %10 from middle half               
        exampleslice_i=ceil(Nsect/2);  %index from middle
        %pick a 'band' of kymographs to analyze
        slices2analyze=round(Ly*linspace(0.5-band/2,0.5+band/2,Nsect));
        ThisKymo_XT.kymopositions=slices2analyze; 
        ThisKymo_XT.cor.prf_tcor=[];  
        ThisKymo_XT.cor.envelope_abs=[];
        
        for ni=1:Nsect  %analyse the kymographs; collect the results                
            slice_i=slices2analyze(ni);
            %1 load, prepare
            resliceim=double(imread(strcat(initval.maindirname,initval.experimentseriesname,'_cln_rs\',reslicefilname),'Index',slice_i));   
            if strcmp(whichaxis, 'y'), resliceim=resliceim'; end            
            resliceim=Clean_Kymo(resliceim);
            
            %2  correlation analysis on Kymograph: peaks, envelope
            if rp==4
                dum=1;
            end
            [cormap_xt_i,t_axis_i,prf_tcor_i,prf_props]=XT_AutoCor(resliceim); 
            
            %Lw=Count_wavelength(resliceim);
            %prf_props.pk2pos=Lw;
            
            
            ThisKymo_XT=Pass_Values(ThisKymo_XT,ni,prf_props,Nsect,'AddValues',fr2min);   
            ThisKymo_XT.cor.prf_tcor=[ThisKymo_XT.cor.prf_tcor ; prf_tcor_i];
            ThisKymo_XT.cor.t_axis=fr2min*t_axis_i; 
            if 0
                [xt_upper,~] = envelope(abs(prf_tcor_i),2,'peak');
                ThisKymo_XT.cor.envelope_abs=[ ThisKymo_XT.cor.envelope_abs; xt_upper]; %add envelope
            end
             %3) boundary count
            % labdadomain=Count_phase_boundaries(resliceim,prf_props);
             labdadomain=NaN;
             ThisKymo_XT.fourier.xtdomains_band(ni)=labdadomain;
                         
             %4) example kymograph & data           
            if ni==exampleslice_i  %set example
                ThisKymo_XT.cor.example.plotim=resliceim;
                ThisKymo_XT.cor.example.cormap_xt=cormap_xt_i;                                               
                        end

        end
        
        %collect averages
        ThisKymo_XT=Pass_Values(ThisKymo_XT,ni,prf_props,Nsect,'AverageValues',fr2min); 
        
        ThisKymo_XT.fourier.xtdomains_av=...
            nanmean(ThisKymo_XT.fourier.xtdomains_band);

  function ThisMov_XY=Get_XY_Movie_properties(filname,initval);
  %Analyze properties of this movie per frame, and averaged
        ThisMov_XY=struct('cor',[]) ;     
        pix2mu=initval.pix2nm/1000;               
        moviepth=strcat(initval.maindirname,initval.experimentseriesname,'_cln\',filname);
        info = imfinfo(moviepth);    
        [ff,~]=size(info);                              
        frames2analyze=min([initval.reps_permovie ff]);
        example_i=ceil(frames2analyze/2);  %index from middle                
        

        for mv=1:frames2analyze %analyze first ten frames (or less)               
            im=double(imread(moviepth,'Index',mv));      
            %Basic processing: 
            [cormap_i, r_axis_i,radialprofile_i,radialprops_i]=RadialAutoCor(im);  
            ThisMov_XY=Pass_Values(ThisMov_XY,mv,radialprops_i,frames2analyze,'AddValues',pix2mu);
            ThisMov_XY.cor.radialprofiles=[ThisMov_XY.cor.radialprofiles; radialprofile_i];          
            %add envelope    
            if 0
                [rpf_upper,~] = envelope(abs(radialprofile_i),2,'peak');
                ThisMov_XY.cor.envelope_abs=[ThisMov_XY.cor.envelope_abs; rpf_upper];  
            end
            ThisMov_XY.cor.r_axis=r_axis_i;
            if mv==example_i
                ThisMov_XY.cor.example.cormap=cormap_i;
                ThisMov_XY.cor.example.plotim=im;
            end
        end
        % movie average: averages of first 10 frames
        ThisMov_XY=Pass_Values(ThisMov_XY,mv,radialprops_i,frames2analyze,'AverageValues',pix2mu);
       
     
        function resliceim_cln=Clean_Kymo(resliceim_o);
            [xx,tt]=size(resliceim_o);
            %equalize x-axis
            mn_prf_x=mean(resliceim_o');
            resliceim=resliceim_o-repmat(mn_prf_x',1,tt);
            %flatten along t-axis
            md_prf_t=median(resliceim);
            axz=1:tt;
            lo=median(md_prf_t)-2*std(md_prf_t);
            hi=median(md_prf_t)+2*std(md_prf_t);
            fitsel=find((lo<md_prf_t)&(md_prf_t<hi));
            prms=polyfit(axz(fitsel),md_prf_t(fitsel),2);
            prfit=prms(1)*axz.^2+prms(2)*axz+prms(3);
            

            resliceim_cln=resliceim-repmat(prfit,xx,1);
            
            if 0  %devlop/test
                figure(481);
                subplot(2,2,1); pcolor(resliceim_o); shading flat; colormap hot;
                subplot(2,2,2); pcolor(resliceim); shading flat; colormap hot;
                subplot(2,2,3); plot(axz,md_prf_t); hold on;
                                plot(axz,prfit,'r-');
                subplot(2,2,4); pcolor(resliceim_cln); shading flat; colormap hot;
                dum=1;
                close(gcf);
            end
            
            function  labdadomain=Count_phase_boundaries(kymo,prf_props);
             %Count number of phase boundaries
%              prf_props
%                 mn1val: -0.4385
%                 mn1pos: 3
%                 pk2val: 0.2285
%                 pk2pos: 5.8129
             kymo=kymo';
             [tt,xx]=size(kymo);
             cpk=prf_props.pk2pos;
             %there should be a reasonable periodicity:
             okpeak=(~isnan(cpk)&cpk>1);   
             if okpeak  
                %get fourier curves per x-position
                mpk=(round(tt/prf_props.mn1pos));      %main common fourierpeak        
                kymo=kymo-(repmat(mean(kymo),tt,1));  %minus average per position
               fim=((fft(kymo)));
               fftreal=abs(fim);
               
               %find fourier peak near comon maximum (obtaine dvia
               %correlation)
               [~,idxes]=max(fftreal(mpk-1:mpk+1,:));  %near maximum
                fftphase=angle(fim);
                phasemain_x=zeros(xx,1);
                for jj=1:xx
                    phasemain_x(jj)=abs(fftphase(idxes(jj)+mpk-1,jj));
                end
                
                %median smooth to keep persistent jumps
                phasemain_x=MedSmooth(phasemain_x,4,'Median');
                
                %count jumps via simple treshold
                diffphase=abs(diff(phasemain_x));
                phasejumptres=pi/3;
                sel=find((diffphase>phasejumptres)&...
                          (diffphase<2));   %avoid periodic phase jumps
                Nb=length(sel);
                
                if 0  %demo plotting
                close all;
                subplot(4,1,1); pcolor(kymo); shading flat; colormap bone;
                title('x-t kymograph')
                xlabel('position, pixels');
                ylabel('time, frames');
                subplot(4,1,2); pcolor(fftreal); shading flat; colormap bone;
                title('x-omega power sepctrum')
                xlabel('position, pixels');
                ylabel('frequency, 1/pixels');
                subplot(4,1,3); plot(phasemain_x); axis tight;
                title('phase per x-position')
                xlabel('position, pixels');
                ylabel('phase');
                subplot(4,1,4); 
                    plot(diffphase); hold on;
                    plot(0*diffphase+phasejumptres,'r-');
                
                axis tight;
                title('difference')
                xlabel('position, pixels');
                ylabel('phase');
                text(xx/2,1,[num2str(Nb) 'jumps']);
                [~]=ginput(1); pause(0.5);
                end
                  else
                      Nb=NaN;
             end
             labdadomain=xx/Nb;
             
            function ThisMov_XY=Pass_Values(ThisMov_XY,mv,radialprops_i,frames2analyze,whattodo,axcalib);
            switch whattodo
                case 'AddValues'
                if mv==1;
                    ThisMov_XY.cor.pk2val_mv=zeros(1,frames2analyze);
                    ThisMov_XY.cor.pk2pos_mv=zeros(1,frames2analyze);
                    ThisMov_XY.cor.mn1val_mv=zeros(1,frames2analyze);
                    ThisMov_XY.cor.mn1pos_mv=zeros(1,frames2analyze);        
                    ThisMov_XY.cor.radialprofiles=[];
                    ThisMov_XY.cor.envelope_abs=[];
                end
                    ThisMov_XY.cor.pk2val_mv(mv)=radialprops_i.pk2val;
                    ThisMov_XY.cor.pk2pos_mv(mv)=axcalib*radialprops_i.pk2pos;
                    ThisMov_XY.cor.mn1val_mv(mv)=radialprops_i.mn1val;
                    ThisMov_XY.cor.mn1pos_mv(mv)=axcalib*radialprops_i.mn1pos;  
                    
                    ThisMov_XY.cor.decay_val_labda_mv(mv)=radialprops_i.decay_val_labda;
                    ThisMov_XY.cor.decaylength_abs_mv(mv)=axcalib*radialprops_i.decaylength_abs;
                    ThisMov_XY.cor.decaylength_sign_mv(mv)=axcalib*radialprops_i.decaylength_sign;
                    ThisMov_XY.cor.decayperiodicityfactor_mv(mv)=radialprops_i.decayperiodicityfactor;
                    ThisMov_XY.cor.strengthofmodulation_mv(mv)=radialprops_i.strengthofmodulation;   
                case 'AverageValues'
                    ThisMov_XY.cor.pk2val_av=nanmedian(ThisMov_XY.cor.pk2val_mv);
                    ThisMov_XY.cor.pk2pos_av=nanmedian(ThisMov_XY.cor.pk2pos_mv);
                    ThisMov_XY.cor.mn1val_av=nanmedian(ThisMov_XY.cor.mn1val_mv);
                    ThisMov_XY.cor.mn1pos_av=nanmedian(ThisMov_XY.cor.mn1pos_mv);
                    
                    ThisMov_XY.cor.decay_val_labda_av=nanmedian(ThisMov_XY.cor.decay_val_labda_mv);
                    ThisMov_XY.cor.decaylength_abs_av=nanmedian(ThisMov_XY.cor.decaylength_abs_mv);
                    ThisMov_XY.cor.decaylength_sign_av=nanmedian(ThisMov_XY.cor.decaylength_sign_mv);
                    ThisMov_XY.cor.decayperiodicityfactor_av=nanmedian(ThisMov_XY.cor.decayperiodicityfactor_mv);
                    ThisMov_XY.cor.strengthofmodulation_av=nanmedian(ThisMov_XY.cor.strengthofmodulation_mv);   
                                       
                    ThisMov_XY.cor=orderfields(ThisMov_XY.cor);
            end
             
                function [matrix_out,regimematrix_mask]=QuickMatrix_it(testcolumn,chanout,EDratio_out,hgt_out,heights,EDratios,matrixoptions);
        LR=length(EDratios); LH=length(heights); 
        regimematrix_mask=ones(LR,LH);
        for hh=1:LH
            thisH=heights(hh);
            for rr=1:LR
                thisR=EDratios(rr);
                sel=find((hgt_out==thisH)&(EDratio_out==thisR)&(strcmp(chanout,'E')));
                sel=find((hgt_out==thisH)&(EDratio_out==thisR));
                if ~isempty(sel),matrix_out(rr,hh)=nanmedian(testcolumn(sel)); end 
            end
        end
        if matrixoptions.flattenregimes
            %1 'traveling waves'
            %2 'amoeba waves'
            %3 'XXL waves'
            %4 'Quasiparticles'
            %5 'Low activity'
            
            regimematrix_mask=6-flipud([...
                [4 4 4 4 1 1];...
                [2 4 4 1 1 1];...
                [2 2 2 3 3 1];...
                [2 2 3 3 3 1];...
                [5 5 5 5 5 5];...
                        ]);
         
            for regime=1:5
                matrix_out(regimematrix_mask==regime)=nanmean(matrix_out(regimematrix_mask==regime));
            end
        end
        matrix_out(isnan(matrix_out))=0;  %for pcolor
        %add empty columns/rows to fool pcolor
        matrix_out=[matrix_out; 0*matrix_out(end,:)]; %extra row
        matrix_out=[matrix_out  0*matrix_out(:,end)]; %extra col
        
        regimematrix_mask=[regimematrix_mask; 0*regimematrix_mask(end,:)+10]; %extra row
        regimematrix_mask=[regimematrix_mask  0*regimematrix_mask(:,end)+10]; %extra col
        
        
        