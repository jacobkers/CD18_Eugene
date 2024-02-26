function Z_info=JK_diagnosis_GetVerticalInfo(Chromosome,chro_stack_raw,chro_stack_decon,chro_pic,initval,cellno,Z_info,cellidx,Ncells); 
       %This function assesses the vertical structure of the domut, by
       %looking for vertical positions of the intensity maximum per contour
       %point. JK2019, for rebuttal 2
       CX=Chromosome.CartesianContourMax_X; 
       CY=Chromosome.CartesianContourMax_Y;
       
       Zprofiles_raw=Get_Zprofiles(chro_stack_raw,CX,CY);
       Zprofiles_decon=Get_Zprofiles(chro_stack_decon,CX,CY);
          
       %maxima
       [Z_MaxPlaneVal_raw,Z_MaxPlaneNo_raw]=max(Zprofiles_raw);       
       [Z_MaxPlaneVal_decon,Z_MaxPlaneNo_decon]=max(Zprofiles_decon);
       
       %FWHM
       FWHM_raw=GetFWHM(Zprofiles_raw);
       FWHM_decon=GetFWHM(Zprofiles_decon);
       
       %local focus error       
       [~,MainPlane]=GetWorkpicFromStack(chro_stack_decon,'FocalPlane');
       FocalPlaneVal_decon=Zprofiles_decon(MainPlane,:);
       meancounts=nanmean(FocalPlaneVal_decon);
       ZLooploss=nanmean(100*abs(Z_MaxPlaneVal_decon-FocalPlaneVal_decon)./...
           meancounts);
       
       Z_info.Plane_av_raw(cellidx)=nanmean(Z_MaxPlaneNo_raw);
       Z_info.Plane_std_raw(cellidx)=nanstd(Z_MaxPlaneNo_raw);
       Z_info.Plane_av_decon(cellidx)=nanmean(Z_MaxPlaneNo_decon);
       Z_info.Plane_std_decon(cellidx)=nanstd(Z_MaxPlaneNo_decon);
       Z_info.FWHM_raw(cellidx)=FWHM_raw;
       Z_info.FWHM_decon(cellidx)=FWHM_decon;
       Z_info.Outofplanelooperrorperc(cellidx)=ZLooploss;
       dum=1;
       
%% plot menu 

set(figure(1),'Visible', 'off');
%figure(1)
%deconvolved

       subplot(2,2,1);
                pcolor(chro_pic); shading flat, colormap hot; hold on;
                plot(CX,CY, 'b-','LineWidth',2);  
                 title('deconvolved');
       subplot(2,2,3); 
           plot(Zprofiles_decon); hold on;
           title('z-profiles')
           xlabel('z-plane no.');
           ylabel('pixel intensity');
       subplot(2,2,4);
           plot(Z_MaxPlaneNo_decon,'ko-');
           title('position')
           xlabel('contour index');
           ylabel('z-plane no max ');
           axis tight
       
           
       outname=strcat('Cell', num2str(cellno,'% 3.0f'));
       outpath='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Paper2018\Rebuttal 2\On 3D aspect\';
       %save figures
       %saveas(gcf,strcat(outpath,'Z_pics per cell\',outname,'.jpg')); 
%        saveas(gcf,strcat(outpath,'Z_pics per cell\',outname,'.svg')); 
       save(strcat(outpath,'Z_pics per cell\',outname,'.mat'),...
           'chro_pic','CX','CY','Zprofiles_decon','Z_MaxPlaneNo_decon');
       pause(0.01);
       close (gcf);
       if (cellidx==Ncells) |(mod(cellidx,10)==0)  %summary
        save(strcat(outpath,'\Z_info.mat'),'Z_info');
       end
       
 function Zprofiles=Get_Zprofiles(stack,CX,CY);       
       [rr,cc,dd]=size(stack); 
       LC=length(CX);
       for ii=1:dd
           plne=stack(:,:,ii); 
           sumval=sum(plne(:));
           plne_smz=JKD2_IM_smoothJK(plne,2);
           stack(:,:,ii)=plne_smz*sumval/sum(plne_smz(:));
       end
              
        %for each contour point, get the profile
        Z_MaxPlaneNo=zeros(LC,1);
        Z_MaxPlaneVal=zeros(LC,1);
        Zprofiles=zeros(dd,LC);
       for ii=1:LC
           prf=squeeze(stack(round(CY(ii)),round(CX(ii)),:)); %simple vertical profile               
           Zprofiles(:,ii)=prf;
       end  
       
    function FWHM=GetFWHM(Zprofs)
        [~,pp]=size(Zprofs);
        Fii=zeros(pp,1);
        for ii=1:pp
            Zprof=Zprofs(:,ii);
            Zprof=Zprof-min(Zprof);
            Fii(ii)=length(find(Zprof>=max(Zprof)/2));
        end
        FWHM=nanmedian(Fii);
        
       
       
       