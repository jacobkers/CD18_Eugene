function Chromosome=ChromosomeRadialSampling(im,QI, sho,initval);
%This function : first performs a COM centering routine;
%Then evaluates the final radial pattern by expansion-fitting each with its own
%radial average. The angle-dependent expansion factor tells us how large
%elliptic deviations are.
% Jacob Kerssemakers 2016 

close all
 
[rr,cc]=size(im);

%1) Get COM center of pattern (first estimate)
%QI.maxradius=rr/2;
[xm,ym]=TrackXY_by_COM_2Dmoment(im);
%2)  Get the radial profiles (running from 0 to 2pi), build a contour line
%and use the coordinates of this line to re-find the COM center
for tt=1:3 
    [rr,cc]=size(im);
     [XX,YY]=meshgrid(1:cc,1:rr);
     Xsamplinggrid=QI.X0samplinggrid+xm;
     Ysamplinggrid=QI.Y0samplinggrid+ym;
     allprofiles=(interp2(XX,YY,im,Xsamplinggrid,Ysamplinggrid,'Linear',0));       
     angles=mod(QI.angles,2*pi);
     [ChromosomeAnnularAxis,ix]=sort(angles);
     Chromosome.RadialMap=allprofiles(ix,:); 
     Chromosome.AnnularAxis=ChromosomeAnnularAxis;
     Chromosome.RadialAxis=QI.radii;
      [aa,~]=size(Chromosome.RadialMap);
     Chromosome.MappingIndex=(1:aa);
 
    %3) Add Cell edge detection--------------------------------
     Chromosome=Get_ChromosomeEdges(Chromosome,QI);       
    %4) Get Cartesian Coordiantes
     angles=Chromosome.AnnularAxis;
     polarmax=Chromosome.PolarContourMax;
     polaredge=Chromosome.PolarContourEdge; 
     CX=xm+polarmax.*cos(angles);
     CY=ym+polarmax.*sin(angles);
            
     %update xm, ym
    [xm,ym,theta,ecc]=JKD2_XY_calculate2Dmomentpoints([CX CY],0); 

end
Chromosome.xCOM=xm;
Chromosome.yCOM=ym; 
Chromosome.CartesianContourMax_X=CX;
Chromosome.CartesianContourMax_Y=CY; 
Chromosome.CartesianContourEdge_X=Chromosome.xCOM+...
                                  polaredge.*cos(angles);
Chromosome.CartesianContourEdge_Y=Chromosome.yCOM+...
                                  polaredge.*sin(angles);                             


 if sho
            figure(4);
            MappingIndex=(1:aa);
            VI=Chromosome.ValidIdx;
            cr=QI.radialoversampling
            pcolor(Chromosome.RadialMap); colormap bone; shading flat; hold on;
            plot(cr*Chromosome.PolarContourMax(VI), MappingIndex(VI), 'ro', 'Linewidth', 2);
            plot(cr*Chromosome.PolarContourEdge(VI), MappingIndex(VI), 'bo', 'Linewidth', 1);
            plot(cr*Chromosome.PolarContourRim(VI), MappingIndex(VI), 'bo', 'Linewidth', 1);
            title('radial map of brightfield with bf-edge overlay');
            xlabel('radial sampling steps');
            ylabel('angular sampling steps');
            [~]=ginput(1);
        end                             
                              
function Chromosome=Get_ChromosomeEdges(Chromosome,QI);
    %Analyze radial map to find edges etc.      
    RadialMap_ori=Chromosome.RadialMap;    
    [rr,cc]=size(RadialMap_ori);       
    RadialMap=RadialMap_ori;
    %1) Get some general properties using a 'softened' version
        DonutMask=repmat(hann(cc)',rr,1);
        MaskedRadialMap=RadialMap_ori.*DonutMask;
        Chro_Stds=nanstd(MaskedRadialMap'); 
        Chro_Std=nanmedian(std(Chro_Stds));
        GapTreshold=0.1*(4*Chro_Std);
        gapidx=find(Chro_Stds<GapTreshold);  %detect gaps
        nongapidx=find(Chro_Stds>=GapTreshold);  %detect gaps
        Gaps=1.0*Chro_Stds<(GapTreshold)';
        
    if 0
        pcolor(RadialMap); colormap hot; shading flat;
        [~]=ginput(1); close(gcf);
    end
    
    %% 2 Find features of typical radial profile
    for ii=1:rr
        prf=RadialMap(ii,:);
        gap=Gaps(ii);
        Profile=Get_ProfileFeatures(prf,gap);
        FWHM(ii)=Profile.FWHM;
        Content(ii)=Profile.Content; 
        MaxC(ii)=Profile.Max; 
        OuterEdge(ii)=Profile.OuterEdge;
        MinC(ii)=Profile.Min;
    end
       
    %padding gaps with neighbouring points OPTIONAL
    if 1
        MaxC=interp1(nongapidx,MaxC(nongapidx), (1:rr), 'linear');
        OuterEdge=interp1(nongapidx,OuterEdge(nongapidx), (1:rr), 'linear');
        MinC=interp1(nongapidx,MinC(nongapidx), (1:rr), 'linear');
        MaxC=JKD1_PRF_smooth(MaxC',4);  %smooth
        OuterEdge=JKD1_PRF_smooth(OuterEdge',4);  %smooth
        MinC=JKD1_PRF_smooth(MinC',4);  %smooth
    end
    
    Chromosome.PolarContourMax=MaxC/QI.radialoversampling;
    Chromosome.PolarContourEdge=OuterEdge/QI.radialoversampling;;
    Chromosome.PolarContourRim=MinC/QI.radialoversampling;;
    Chromosome.PolarContourContent=Content;
    Chromosome.PolarContourFWHM=FWHM;
    Chromosome.Gaps=Gaps;
    Chromosome.GapIdx=gapidx;
    Chromosome.ValidIdx=nongapidx;

 
    if 0
    figure;    
        plot(Chro_Stds); hold on;
        plot(gapidx,Chro_Stds(gapidx),'rx');   
        [~]=ginput(1); close(gcf);
    end
    

function Profile=Get_ProfileFeatures(prf,gap);
%Get typical radial features of profile
prf=prf-min(prf);
cc=length(prf); 
if ~gap   
    %1) find first profile max from center pos
    prf_han=prf.*(hann(cc)');
    [~,MaskCenter]=max(prf_han);
    MaskWidth=50;
    prf_smz=JKD1_PRF_smooth(prf',20);     
    maxidx=Get_GuidedMax(prf_smz',MaskCenter,MaskWidth);
    peakval=prf(maxidx);
    if 0
        figure;    
        plot(prf,'b-'); hold on;
        plot(prf_smz,'r-');   
        plot(maxidx,peakval,'ko');  
        [~]=ginput(1); close(gcf);
    end
    
    outprf=prf(maxidx:end);  %profile on inside 
    inprf=prf(1:maxidx-1);   %profile on outside
      
    %dark rim
    [~,ix]=min(outprf);
    if isempty(maxidx)
        dum=1;
    end

    
    outminidx=ix+maxidx-1;    
    %first above zero
    [~,ix2]=min(abs(prf(maxidx:outminidx)));
    edgeidx=ix2+maxidx-1;
    
    %content
    sel=find(prf>0);
    cnt=nansum(prf);
    
    %width
    
    
    sel2=find((prf>0.5*peakval)&(~isnan(prf)));
    FWHM=length(sel2);
    
    Profile.FWHM=FWHM;
    Profile.Content=cnt; 
    Profile.Max=maxidx;
    Profile.Min=outminidx;       
    Profile.OuterEdge=edgeidx;

else
    Profile.FWHM=NaN;
    Profile.Content=NaN; 
    Profile.Max=NaN;
    Profile.Min=NaN;
    Profile.OuterEdge=NaN;
    
end 
    


function outdata=JKD1_PRF_smooth(data,span)
%running window 'span' smoothing average; 0=no smoothing.
%perioc boundaries
%JacobKers 2013

if span>0
    hs=ceil(span/2);
    le=length(data);
    outdata=zeros(le,1);
    midpt=ceil(le/2);
    data2=[data(midpt:end) ; data ; data(1:midpt)];
    lex=length(data(midpt:end));
    for i=1:le
    idx=i+lex;
    outdata(i)=mean(data2(idx-hs:idx+hs));
    end
else
    outdata=data; %no smoothing
end

function data_expanded=Expand_elliptoidrange(data);
    %Expand for later interpolation
    [rr,~]=size(data);
    halfidx=ceil(rr/2);
    half1=data(1:halfidx-1,:);
    half2=data(halfidx:end,:);
    half1(:,1)=half1(:,1)+2*pi;
    half2(:,1)=half2(:,1)-2*pi;
    data_expanded=[half2 ; data ; half1];



           
 
       



