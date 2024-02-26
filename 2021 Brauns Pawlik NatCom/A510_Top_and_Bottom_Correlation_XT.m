function A510_Top_and_Bottom_Correlation_XT
JWJK_A:-------------------------------------------------------------------
%Title: load 2 stacks and correlate them image by image; allows for  timeshift
%
%Summary: This code is used to compare MinED patterns facing each other on the top and the bottom of a sample 

%%
%Input: movies taken at different sample heights are loaded. To ave time,
%these movies were tilted and resliced, so that xt-profiles can be loaded at once. 



%Output: correlation data, in mat, excel and jpg

%Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
%code designed & written Jacob Kerssemakers 2016 
%:JWJK_A-------------------------------------------------------------------
%

codepth='D:\jkerssemakers\_Recent\CD\BN_CD16_Greg\GregCodeExchangePack\MatlabCode';
outpath='D:\jkerssemakers\_Recent\CD\BN_CD16_Greg\20181009 TopBottom\ViaTimeTraces\';
addpath(genpath(codepth));
close all;

showplots=1;
%1 load single image

%% data 
experimentseriesname='movies_test';  %the name used for this series
exampledirname='D:\jkerssemakers\_Data\CD\2016_Greg\20181009_bottomtop_cleaned\';
maindirname='D:\jkerssemakers\_Data\CD\2016_Greg\20181009_bottomtop_cleaned\Resliced\';
heights=[10 20 21 26 30 43 51];
%heights=[20 43];

LH=length(heights);

cordata=zeros(LH,5);
allcorvals=[];
allcorhist=[];
for ii=1:LH;
        hgt=heights(ii)
        top_pth=strcat(maindirname,num2str(hgt),'top_cleaned_resliced_x.tif');
        bot_pth=strcat(maindirname,num2str(hgt),'bottom_cleaned_resliced_x.tif');
        top_ex=strcat(exampledirname,num2str(hgt),'top_cleaned.tif');
        bot_ex=strcat(exampledirname,num2str(hgt),'bottom_cleaned.tif');

        first_top=double(imread(top_pth,'Index',1));
        info = imfinfo(top_pth);    
        [yy,~]=size(info);                                           
        [xx,tt]=size(first_top);  
        
        top_example=double(imread(top_ex,'Index',10));  
        bot_example=double(imread(bot_ex,'Index',10)); 
        
        cormap=zeros(yy,xx);   
        
        crop1st=5;  %remove first frames (sometimes artefacts)
        for yi=1:yy  %for every 'kymo-slab'
            top=double(imread(top_pth,'Index',yi));  
            bot=double(imread(bot_pth,'Index',yi)); 
            for xi=1:xx
                top_xt=top(xi,:);
                bot_xt=bot(xi,:);
                top_xt=Clean_Normalize_it(top_xt(crop1st:end));
                bot_xt=Clean_Normalize_it(bot_xt(crop1st:end));
                corcurve=(real(ifft(fft(top_xt).*conj(fft(bot_xt))))); 
                corcurve_auto=(real(ifft(fft(bot_xt).*conj(fft(bot_xt)))));        
                cormap(yi,xi)=corcurve(1)/corcurve_auto(1); %zero-shift correlation               
%                 plot(corcurve); hold on;
%                 plot(corcurve_auto,'r-'); hold off;
%                 pause(0.3); close(gcf);
            end
        end
       plotcormap=cormap;      
       plotcormap(1:end,1:2)=1; plotcormap(1:end,3:4)=-1;  %scale marker
       cormap4hist=cormap(5:end,:);  %remove inactive strip
       binaxis=-1:0.1:1;
       cormap4histdata=cormap4hist(:);
       corhist=hist(cormap4histdata,binaxis);
       corhistperc=100*corhist/sum(corhist);
       allcorvals=[allcorvals  cormap4hist(:)];
       allcorhist=[allcorhist  corhistperc'];
       
       %some further analysis
       s1=find(binaxis<=-0.6);
       s2=find(abs(binaxis)<0.6);
       s3=find(binaxis>=0.6);
       perc_anticor=sum(corhistperc(s1));
       perc_uncor=sum(corhistperc(s2));
       perc_cor=sum(corhistperc(s3));
       
       cordata(ii,:)=[hgt perc_cor perc_anticor ...
                      perc_cor+perc_anticor perc_uncor]; 
   
       if showplots
           close(gcf);
           
            subplot(2,2,1); pcolor(top_example); colormap hot; shading flat;
            title(['top',num2str(hgt)]);
            subplot(2,2,2); pcolor(bot_example); colormap hot; shading flat;
            title(['bottom',num2str(hgt)]);
            subplot(2,2,3); pcolor(plotcormap); colormap hot; shading flat;
            title('Correlation');
            subplot(2,2,4); 
                bar(binaxis,corhist,'b');
                title('local top-bottom correlation');
                xlabel('correlation value, a.u.');
                ylabel('area counts, a.u.');
            pause(0.02);
            matrixname=strcat('Height',num2str(hgt));
            figoutname=strcat(outpath,matrixname,'_overview','.jpg');
            dataoutname=figoutname(1:end-4);
            saveas(gcf,figoutname, 'jpg');  
            %save raw data
            save(dataoutname,'cormap4hist','plotcormap','cormap4histdata','binaxis');
            dum=1;
        end
end
figure(2);
plot(cordata(:,1),cordata(:,2:end-1),'o-', 'Linewidth',2);
xlabel('height, microns');
ylabel('area coverage, %');
legend('correlated', 'anticorrelated', 'total', 'uncorrelated')

% subplot(2,1,2);
%     [rr,cc]=size(allcorvals);
%     xax=repmat(heights,rr,1);
%     scatax=randn(rr,length(heights))
%     plot(xax+scatax,allcorvals, 'x', 'MarkerSize', 1)
%     xlabel('height, microns');
%     ylabel('local correlation value, a.u.')

dum=1;

figoutname=strcat(outpath,'A510_Top-bottom_Correlation_heightdependence.jpg');
            saveas(gcf,figoutname, 'jpg');  

%%organize export
%1) Excel file

xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),{'heights'},'Histograms','B1');
xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),heights,'Histograms','B2');
xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),{'Correlation value'},'Histograms','A2');

xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),[binaxis' allcorhist],'Histograms','A3');       

ColNames=[{'height'} {'correlated area %'} {'anticorrelated area %'}...
           {'all (anti)correlated area %'} {'uncorrelated area %'}];
xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),ColNames,'Categorized','A1');
xlswrite(strcat(outpath,'A510_Top-bottom_Correlation_results.xlsx'),cordata,'Categorized','A2');


   
 function prf=Clean_Normalize_it(prf_in);
     %flatten, smooth and normalize a curve
     LP=length(prf_in);
     axz=1:LP;
     prf=prf_in-mean(prf_in); 
     params=polyfit(axz,prf,1);
     prfit=axz*params(1)+params(2);
     prf=prf-prfit;
     prf=smooth(prf,2);
     prf=(prf-mean(prf))/(4*std(prf));
%       plot(axz,prf,'r'); hold off;
%      dum=1;
%         if showplots
%             subplot(2,2,4); plot(RelCor);             
%             pause(0.5);    
%             matrixname=strcat('Height',num2str(hgt),'DeltaT',num2str(shft));
%             figoutname=strcat(outpath,matrixname,'_overview_','.jpg');
%             saveas(gcf,figoutname, 'jpg');
%             close(gcf);
%         end

% 
% outdata=[timeshift allcorrelations];   
%       



