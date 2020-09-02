function kymo=kym_build_kymo_from_movie(pth,expinfo)
%Build a kymograph from a movie stack or sequence; 
%use start  and endpoints, width and drift
%JacobKers2019--------------------------------------

    codepth=pwd;   
    X=expinfo.endpoints_xy(:,1);
    Y=expinfo.endpoints_xy(:,2);
    driftx=expinfo.driftxy(1); 
    drifty=expinfo.driftxy(2);
    if isfield(expinfo, 'pad_it');
        pad_it=expinfo.pad_it;
    else
        pad_it=0;
    end
    if pad_it
         padL=10;  Y(2)=Y(2)+2*padL;
    end
    hf=ceil(expinfo.kymowidth/2);
        
    cd(pth); fs=dir('*.tif*'); cd(codepth);
    %1 check directory, build filename
     [filno,~]=size(fs);
     if filno>1, option='single_tifs'; else option='stack';end
     filname1=fs(1).name;
     firstplane=double(imread(strcat(pth,'\',filname1),'Index',1));
     info = imfinfo(strcat(pth,'\',filname1));  
     filnametemplate=filname1(1:end-6);
     [rr,cc]=size(firstplane);     
     switch option
         case 'stack', [FramesNo,~]=size(info);
         case 'single_tifs',[FramesNo,~]=size(fs);
     end
      %2) work all images
     for idx=1:FramesNo    
         if mod(idx,250)==0, disp(num2str(idx)),end
        switch option
         case 'single_tifs' 
             filname_i=[filnametemplate, num2str(idx),'.tiff'];
             %180514_203637-1
             pic=double(imread(strcat(pth,'\',filname_i)));
         case 'stack'      
             pic=double(imread(strcat(pth,'\',filname1),'Index',idx));
        end      
        
        %add optional padding here.....
        if pad_it
            %[outpic,thr,SR_ratio]=matrix_find_treshold_MD2020(pic,0);
            noiz=std(pic(:));
            av=median(pic(:));
            padblock=av+noiz*randn(padL,cc);
            pic=[padblock ; pic; padblock];
            [rr,cc]=size(pic);
            %get noise       
        end
        
        %drift-corrected sampling
        Xi=X+driftx*idx/FramesNo;
        Yi=Y+drifty*idx/FramesNo;
        
        prf=xy_get_wide_intensityprofile(pic,Xi,Yi,hf);        
        if idx==1
            kymo=zeros(FramesNo,length(prf));
            [~,Wk]=size(kymo);
        end     
        
        prf=prf_remove_sloped_background(prf);          
        lp=min(Wk,length(prf));% avoid rounding effects       
        kymo(idx,1:lp)=prf(1:lp); 
     end
     dum=1;
     