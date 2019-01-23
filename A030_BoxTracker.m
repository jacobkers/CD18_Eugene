function A030_BoxTracker
close all;
figure(1);

loaddatatype='Sim';
loaddatatype='ASCII';
initval.boxhalfwidth=2;
initval.tracklookahead=20;
initval.smoothlookahead=10;
clickit=1;

switch loaddatatype
    case 'Sim'  %out of function
        trackmap=Simulate_it;        
    case 'ASCII'
        expname='2019_01_22 non_interactive type';
        % [no xL tL1 tL2 xR tR1 tR2]              
        roistartstop=[ [9 68.155 2940 1E6 68.343 300 1E6];...
                       [31 65.25 432 1E6 67.726 1450 1E6];...
            ];
        datapath='D:\jkerssemakers\_Data\CD\2018_Eugene\';
        exppath=[datapath,'\',expname,'\'];
        outpath=strcat(datapath, 'matlabresults\',expname,'\');
        textfile='Kymograph_DNA.txt';
        if ~isdir(outpath), mkdir(outpath); end
end


[Lroi,~]=size(roistartstop);
for roii=1:Lroi
    close all;
    roino=roistartstop(roii,1);
    Exp=strcat('ROI',num2str(roino));
    SaveName=char(strcat(outpath, Exp));
    datainpath=strcat(exppath,'M', num2str(roino),'\kymo_ImageJ\');       
    source=[ datainpath, textfile];
    trackmap=dlmread(source);
    [ff,cc]=size(trackmap);
    savit=1;
 
    pcolor(trackmap); colormap hot, shading flat; hold on;
    if clickit
    clickpoints=ginput(4);  %click [no xL tL1 tL2 xR tR1 tR2]% 
    thisroistartstop=[
            [roino  clickpoints(1,1) clickpoints(1,2) clickpoints(2,2)...
                    clickpoints(3,1) clickpoints(3,2) clickpoints(4,2)] 
            ];        
    else
        thisroistartstop=roistartstop(roii,:);
    end
    %clean it
    if thisroistartstop(4)>ff-initval.tracklookahead,
        thisroistartstop(4)=ff-initval.tracklookahead;
    end
    if thisroistartstop(7)>ff-initval.tracklookahead,
        thisroistartstop(7)=ff-initval.tracklookahead;
    end   
    
    colormap(bone);
    plot(thisroistartstop(2),thisroistartstop(3),'ro-', 'MarkerFaceColor', 'r');
    plot(thisroistartstop(5),thisroistartstop(6),'bo-', 'MarkerFaceColor', 'b'); 
    plot([5 cc-5]',[thisroistartstop(4) thisroistartstop(4)]','wo-');
    plot([5 cc-5]',[thisroistartstop(7) thisroistartstop(7)]','wo-');
    
    pause(0.5);
    close(gcf);
 
    initval.tetherlevel=Get_tetherlevel(trackmap);
    for jj=1:2 %initialize loop info
        xj=thisroistartstop(2+(jj-1)*3); %first x for trace;
        fj=thisroistartstop(3+(jj-1)*3); %first frame;
        [~,Ij,~]=Get_box_intensity(trackmap,xj,fj,initval,0); %first I
        looptraces(jj).Lx=xj;
        looptraces(jj).frame=fj;
        looptraces(jj).LI=Ij;
    end 
    
    for ii=1:ff
        %if mod(300,ii), disp(num2str(ii));end;
        looptraces=Analyze_traces(looptraces,trackmap,thisroistartstop,ii,initval);
    end


    subplot(1,2,1);
    pcolor(trackmap); colormap bone, shading flat; hold on;
    colormap(bone);
    plot(looptraces(1).Lx+0.5,looptraces(1).frame+0.5,'r-', 'MarkerFaceColor', 'w');
    plot(looptraces(2).Lx+0.5,looptraces(2).frame+0.5,'b-', 'MarkerFaceColor', 'w');
    ylabel('frameno, a.u.');
    xlabel('position, a.u.');

    subplot(1,2,2);
    plot(looptraces(1).frame,looptraces(1).LI,'r-', 'MarkerFaceColor', 'w'); hold on;
    plot(looptraces(2).frame,looptraces(2).LI,'b-', 'MarkerFaceColor', 'w'); hold on;
    xlabel('frameno, a.u.');
    ylabel('loop intensity, % of total');
    legend('loop 1', 'loop 2');


    if savit
    save(strcat(SaveName, '_pairtrackresults.mat'),... 
                                  'looptraces');

    saveas(gcf, strcat(SaveName, '_pairtrackresults.svg'));
    saveas(gcf, strcat(SaveName, '_pairtrackresults.jpg'));
    end
end

function looptraces=Analyze_traces(looptraces,trackmap,thisroistartstop,ii,initval);
    
    %[no xL tL1 tL2 xR tR1 tR2]% 
    hf=initval.boxhalfwidth;    
    for iL=1:2 %for each loop
        loopexist=(ii>=thisroistartstop(3+(iL-1)*3)&(ii<thisroistartstop(4+(iL-1)*3)));
        if loopexist  %update the loop position
            i_cur=length(looptraces(iL).Lx); %next loop index            
            x_cur=looptraces(iL).Lx(i_cur);
            fr_nxt=looptraces(iL).frame(i_cur)+1;
            [~,~,box]=Get_box_intensity(trackmap,x_cur,fr_nxt,initval,initval.tracklookahead); %get next box
            
            %new position
            bx=box-min(box);
            mbox=GaussMask(box,0.5);
            [com,comc]=Get_1DCOM(mbox);
            x_nxt=x_cur+comc;   %next position
            
            %for all other loops: check minimum distance
            
%             if 0   %add auto-repulsion               
%                 if (abs_Rx-abs_Lx)<2*hf
%                     repulsion=(2*hf-(abs_Rx-abs_Lx))/2;
%                     abs_Rx=abs_Rx+repulsion;
%                     abs_Lx=abs_Lx-repulsion;
%                 end
%             end            
            [~,Iperc,~]=Get_box_intensity(trackmap,x_nxt,fr_nxt,initval,initval.smoothlookahead);           
            
            %update trace
            looptraces(iL).frame(i_cur+1)=fr_nxt;
            looptraces(iL).Lx(i_cur+1)=x_nxt;
            looptraces(iL).LI(i_cur+1)=Iperc;  

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

         
  function [Iperc,Iperc_excess,box]=Get_box_intensity(trackmap,xx,fr,initval,lookahead);
        [~,cc]=size(trackmap);
        fr=round(fr);
        hf=initval.boxhalfwidth;      
        lo=round(max([1 xx-hf])); 
        hi=round(min([cc xx+hf]));        
        if lookahead>0
            squbox=trackmap(fr:fr+lookahead,lo:hi);
            maskline=fliplr(linspace(0.1,1,lookahead+1));
            maskline=maskline/sum(maskline);
            squmask=repmat(maskline',1,hi-lo+1);
            box=sum(squbox.*squmask);
        else
            box=trackmap(fr,lo:hi);
        end
        Irw=sum(box);
        Irw_excess=Irw-initval.tetherlevel*length(box);
        Itotal=sum(trackmap(fr,:));
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
    
 function tetherbackgroundlevel=Get_tetherlevel(trackmap);   
    [rr,cc]=size(trackmap);
     midpart=trackmap(:,round(cc/4):round(3*cc/4));
     tetherbackgroundlevel=median(midpart(:));
     dum=1;
    
  function trackmap=Simulate_it;
        initval.BallParkRadius=1;
        initval.PadCurves=1;
        LT=100;
        LX=50;
        kymo=2*rand(LT,LX);
        Taxis=1:LT;
        Spotpos1=round((LX/3)*Taxis/LT+LX/3+1*rand(1,LT));
        Spotpos2=round(Spotpos1+5+1*randn(1,LT));
        for ii=1:LT
           for jj=-3:3
           val1=kymo(Taxis(ii),Spotpos1(ii)+jj);
           kymo(Taxis(ii),Spotpos1(ii)+jj)=val1+(3-abs(jj))^2;

           val2=kymo(Taxis(ii),Spotpos2(ii)+jj);
           kymo(Taxis(ii),Spotpos2(ii)+jj)=val2+(3-abs(jj))^2;
           end
        end
        trackmap=kymo;  
    