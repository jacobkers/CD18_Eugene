function A030_BoxTracker
close all;

expname='Figure2Pannel'; 
initval=A001_Initialize_BoxTracker(expname);
[~,Lroi]=size(initval.roistartstop);

%%  work trhough all ROIs
for roii=1:Lroi
    close all;
    figure(roii);
    
    %% load a region-of-interest kymograph
    savit=1;
    roino=initval.roistartstop(roii).roino;
    Exp=strcat('ROI',num2str(roino));
    SaveName=char(strcat(initval.expi_outpath, Exp));
    datainpath=strcat(initval.expi_inpath,'M', num2str(roino),'\kymo_ImageJ\');       
    source=[datainpath, initval.kymofile];
    trackmap=dlmread(source);
    %% get general properties, such as tether edges
    [ff,cc]=size(trackmap);
    thisroistartstop=initval.roistartstop(roii);   
    [tetherstart,tetherstop]=Get_tetheredges(trackmap);
    initval.tetherlevel=Get_tetherlevel(trackmap,tetherstart,tetherstop);
    initval.tetherstart=tetherstart;
    initval.tetherstop=tetherstop;       
    Nloops=length(thisroistartstop.startx);
    looptraces=struct('Lx',[]);
    
    %% initialize loop info
    for jj=1:Nloops 
        xj=thisroistartstop.startx(jj); %first x for trace;
        fj=thisroistartstop.pre_t(jj); %first frame;
        %get first box
        boxprops.boxhalfwidth=thisroistartstop.loopanalysishalfwidth(jj);
        boxprops.cutlevel=initval.tetherlevel;
        boxprops.lookahead=0;
        [~,Ij,~]=Get_box_intensity(trackmap,xj,fj,boxprops); %first I
        looptraces(jj).Lx=xj;
        looptraces(jj).frame=fj;
        looptraces(jj).I_mid=Ij;
    end
    
    %% now, analyze the kymograph
    for ii=1:ff
        if (mod(ii,300)==0), disp(num2str(ii));end;
        looptraces=Analyze_traces(looptraces,trackmap,thisroistartstop,ii,initval);
    end

    %% Plot menu
    subplot(1,2,1);
    pcolor(trackmap); colormap bone, shading flat; hold on;
    colormap(bone); 
    for jj=1:Nloops 
        plot(looptraces(jj).Lx+0.5,looptraces(jj).frame+0.5,'r-');
        %length()
        plot(looptraces(jj).curvestart+0.5,looptraces(jj).frame+0.5,'y-');
        plot(looptraces(jj).curvestop+0.5,looptraces(jj).frame+0.5,'y-');
        ylabel('frameno, a.u.');
        xlabel('position, a.u.');
    end
    subplot(1,2,2);
    for jj=1:Nloops
        plot(looptraces(jj).frame,looptraces(jj).I_mid,'y-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).I_left,'m-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).I_right,'b-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).checksum,'k-'); hold on;
    end
%         checksum=looptraces(jj).I_loop+...
%                   looptraces(jj).I_left+...
%                   looptraces(jj).I_right;
%         plot(looptraces(jj).frame,checksum,'b-'); hold on;
        legend('loop','left','right');
        xlabel('frameno, a.u.');
        ylabel('loop intensity, % of total');
    

    if savit
    save(strcat(SaveName, '_pairtrackresults.mat'),... 
                                  'looptraces','initval');
    saveas(gcf, strcat(SaveName, '_pairtrackresults.svg'));
    saveas(gcf, strcat(SaveName, '_pairtrackresults.jpg'));
    saveas(gcf, strcat(SaveName, '_pairtrackresults.fig'));
    end
end

function looptraces=Analyze_traces(looptraces,trackmap,thisroistartstop,ii,initval);   
    %[no xL tL1 tL2 xR tR1 tR2]% 
    Nloops=length(thisroistartstop.startx);   
    for iL=1:Nloops %for each loop
        crp=max([initval.smoothlookahead initval.tracklookahead]);
        loopstandstart=thisroistartstop.pre_t(iL);
        loopwalkstart=thisroistartstop.startt(iL);
        loopwalkstop=thisroistartstop.stopt(iL)-crp;
        %% perform track analysis in predefined sections of loop
        loopalife=(ii>=loopstandstart&(ii<loopwalkstop));
        loopstands=(loopalife&(ii<loopwalkstart));
        loopmoves=(loopalife&(ii>=loopwalkstart));
      
        
        %% analyze the pre-set life slot of the loop   
        if loopalife   
            i_cur=length(looptraces(iL).Lx); %next loop index            
            x_cur=looptraces(iL).Lx(i_cur);
            fr_nxt=looptraces(iL).frame(i_cur)+1;
            boxprops.boxhalfwidth=thisroistartstop.trackhalfwidth(iL);
            boxprops.cutlevel=initval.tetherlevel;
            boxprops.lookahead=initval.tracklookahead;
            [~,~,box]=Get_box_intensity(trackmap,x_cur,fr_nxt,boxprops); %get next box
            box(box<initval.tetherlevel)=initval.tetherlevel;  %padding off-tether
            box=box-min(box);
            %new position
            
            %box=box-min(box);
            mbox=GaussMask(box,0.5);
            [com,comc]=Get_1DCOM(mbox);
            %% analyze the pre-loop time slot (were intensity is absent or low)
            if loopstands
                x_nxt=x_cur;   %next position; do not update
            end
            %% analyze the motion time slot (where the loop is active)
            if loopmoves
                x_nxt=x_cur+comc;   % update next position; no update!
            end

            %% get intensity of newly tracked box; more smoothened
            boxprops.boxhalfwidth=thisroistartstop.loopanalysishalfwidth(iL);
            boxprops.cutlevel=initval.tetherlevel;
            boxprops.lookahead=initval.smoothlookahead;
            [~,~,box_sm,boxlo,boxhi,box_full]=Get_box_intensity(trackmap,x_nxt,fr_nxt,boxprops);           
            [strt,stp,~]=Get_edgeslength(box_sm,'loop');
            loopedge_lft=boxlo-1+strt;
            loopedge_rgt=boxlo-1+stp;
            
            %Get corrected left and right intensities
            %new position from smoothened box
            box_sm(box_sm<initval.tetherlevel)=initval.tetherlevel;  %padding off-tether
            box_sm=box_sm-min(box_sm);
            mbox_sm=GaussMask(box_sm,0.5);
            [com_sm,~]=Get_1DCOM(mbox_sm);           
            mid=com_sm-1 +boxlo;
         
            box_res=box_full;
            box_res(boxlo:boxhi)=box_res(boxlo:boxhi)-box_sm;
            
            %note:sum(box_res)+sum(box_sm)=sum(box_full)
            
            box_left=box_res(1:loopedge_lft-1);
            box_stem=box_res(loopedge_lft:loopedge_rgt);
            box_right=box_res(loopedge_rgt+1:end);
            
            %note: sum(box_res)=sum(box_left)+sum(box_stem)+sum(box_right)
            
            looptraces(iL).frame(i_cur+1)=fr_nxt;
            looptraces(iL).curvestart(i_cur+1)=boxlo-1+strt;
            looptraces(iL).Lx(i_cur+1)=x_nxt;           
            looptraces(iL).curvestop(i_cur+1)=boxlo-1+stp;
            looptraces(iL).I_left(i_cur+1)=100*sum(box_left)/sum(box_full);
            looptraces(iL).I_right(i_cur+1)=100*sum(box_right)/sum(box_full);
            looptraces(iL).I_stem(i_cur+1)=100*sum(box_stem)/sum(box_full);
            looptraces(iL).I_hat(i_cur+1)=100*sum(box_sm)/sum(box_full);;
            looptraces(iL).I_mid(i_cur+1)= looptraces(iL).I_hat(i_cur+1)+...
                                           looptraces(iL).I_stem(i_cur+1);
            looptraces(iL).checksum(i_cur+1)=...
                    looptraces(iL).I_left(i_cur+1)+...
                    looptraces(iL).I_mid(i_cur+1)+...
                    looptraces(iL).I_right(i_cur+1);                           
            if i_cur+1>250
                dum=1;
            end
            %update trace
            

            

             if 0% iL==2
             figure(2);
             plot(box,'b-'); hold on; 
             plot(mbox,'r-'); hold on;  
            stem(com,max(box),'ro');
             hold off;
            pause(0.2);
            [~]=ginput(1);
            end
        end
    end
     
    function [com,comc]=Get_1DCOM(soi);
        %Get one-dimensional center of mass (as offset from center)
        lp=length(soi);
        ax=linspace(-lp/2,lp/2,lp)';
        mn=nanmin(soi);
        if ~isempty(mn),soi2=soi'-mn; else soi2=soi'; end
        %background correction
        comc=sum(soi2.*ax)/sum(soi2);
        com=comc+lp/2+0.5;      %(center of mass in array coordinates)
        com=max([1 com]); com=min([lp-1 com]); %just to be sure
        comc=max([-lp/2 comc]); comc=min([lp/2 comc]); %just to be sure

         
  function [Iperc,Iperc_excess,box,lox,hix,box_full]=Get_box_intensity(trackmap,xx,fr,boxprops);
        [tt,cc]=size(trackmap);
        fr=round(fr);
        hf=boxprops.boxhalfwidth;      
        lox=round(max([1 xx-hf])); 
        hix=round(min([cc xx+hf]));
        lofr=round(max([1 fr])); 
        hifr=round(min([tt fr+boxprops.lookahead]));
        
        if boxprops.lookahead>0
            squbox=trackmap(lofr:hifr,lox:hix);
            squbox_full=trackmap(lofr:hifr,:);
            
            maskline=fliplr(linspace(0.1,1,hifr-lofr+1));
            maskline=maskline/sum(maskline);
            
            squmask=repmat(maskline',1,hix-lox+1);
            squmask_full=repmat(maskline',1,cc);
            
            box=sum(squbox.*squmask);
            box_full=sum(squbox_full.*squmask_full);
        else
            box=trackmap(fr,lox:hix);
            box_full=trackmap(fr,1:cc);
        end
        box_full=box_full-min(box_full);
        
        %to avoid dark background issues
        box(box<boxprops.cutlevel)=boxprops.cutlevel;
             
        Irw=sum(box);
        Irw_excess=Irw-boxprops.cutlevel*length(box);
        
        Itotal=sum(box_full);
        Iperc=100*Irw/Itotal;
        Iperc_excess=100*Irw_excess/Itotal;

function ar=GaussMask(ar,sigma)
    %This function masks an array with a Gaussian window;
    %sigma=1 means edge is on one sigma etc.
    lar=length(ar);
    x0=lar/2;
    sig=sigma*x0;
    ax=linspace(-x0,x0,lar);
    ar=ar.*exp(-(ax/sig).^2);
    
 function tetherbackgroundlevel=Get_tetherlevel(trackmap,tetherstart,tetherstop);   
    [rr,cc]=size(trackmap);
     midpart=trackmap(:,tetherstart+5:tetherstop-5);
     tetherbackgroundlevel=median(midpart(:));
     dum=1;
    
  
        
     function [tetherstart,tetherstop]=Get_tetheredges(trackmap)
        [ff,cc]=size(trackmap);
            %get the average tether position (assume it is fixed)
        all_tetheredges=zeros(ff,2);
        for jj=1:ff  
            prf_res=trackmap(jj,:);
            prf_res=prf_res-min(prf_res);
            [t_strt,t_stp,~]=Get_edgeslength(prf_res,'tether');
            all_tetheredges(jj,:)=[t_strt t_stp];
        end
        tetherstart=nanmedian(all_tetheredges(:,1));
        tetherstop=nanmedian(all_tetheredges(:,2));
          dum=1; 
        

 function [curvestart,curvestop,curveok]=Get_edgeslength(prf_res, type_of_profile);
    %function uses 'shaved off' profile to find start and stop
    tresval=0.4;
    Lp=length(prf_res); axz=1:Lp;
        
    [lo_L,ixL]=min(prf_res(1:ceil(Lp/2)));
    [lo_R,ixR]=min(prf_res(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
    slopefit=polyval(polyfit([ixL ixR],[lo_L lo_R],1),axz);
    
    prf_res_ft=prf_res-slopefit;
    
    switch type_of_profile
        case 'tether' %main level, assuming most of profile is tether 
            main_lev=median(prf_res_ft);
            %subplot(2,2,2);
        case 'loop' %%FWHM level, assuming relatively compact loop 
            main_lev=0.5*max(prf_res_ft);
            %subplot(2,2,4);
    end
    lo=min(prf_res_ft);
    sel=find(prf_res_ft>(lo+tresval*(main_lev-lo)));
    
    if ~isempty(sel);
        curvestart=min(sel);     
        curvestop=max(sel);
        curveok=1;
    else
        curvestart=NaN;     
        curvestop=NaN;
        curveok=0;
    end

    if 0      
        plot(prf_res,'b-'); hold on;
        plot(slopefit,'b-');
        plot(prf_res_ft,'r');       
        plot(0*prf_res+main_lev,'r-');
        title(type_of_profile);
        stem(curvestart,prf_res_ft(curvestart),'ro');
        stem(curvestop,prf_res_ft(curvestop),'ro');
        legend( '\fontsize{6} residu','\fontsize{6} minima fit',...
                '\fontsize{6} corrected','\fontsize{6} mean',...
                '\fontsize{6}edges');
        xlabel('positition, pixel units');
        ylabel('fluorescence intensity, a.u.');
        pause(0.5);        
        [~]=ginput(1);
        hold off;
    end         
          
 