function G00_ProcessFurther_Get_LengthOfContentAroundLocation(Axis,Content,PercentualShiftFromCfpPos,PercOfContent)
%This function calculates the length of a pre-set percentage of genomic
%content around a relative position 
blowup=10;  %increase resolution, equidistant points
PosPerc=100-initval.Rfppos+initval.Cfppos-PercentualShiftFromCfpPos;
BlowUpAxis=linspace(Axis(1), Axis(end),10*length(Axis));
BlowUpContent=interp1(Axis,Content,BlowUpAxis);

%BlowUpContent=BlowUpContent/sum(BlowUpContent)*100;
%Normalize (after putting points equidistant)

[~,StartIdx]=min(abs(BlowUpAxis-PosPerc))

AxLeft=fliplr(PosPerc-BlowUpAxis(1:StartIdx));  %split axis
AxRight=BlowUpAxis(StartIdx:end)-PosPerc;

ContentLeft=fliplr(BlowUpContent(1:StartIdx));  %content per pos
ContentRight=BlowUpContent(StartIdx:end);

%Scale Contents
CumContentLeft=cumsum(ContentLeft-ContentLeft(1));
CumContentLeft=CumContentLeft/CumContentLeft(end)*AxLeft(end);
CumContentRight=cumsum(ContentRight-ContentRight(1));
CumContentRight=CumContentRight/CumContentRight(end)*AxRight(end);


[~,LeftIdx]=min(abs(CumContentLeft-PercOfContent/2));
[~,RightIdx]=min(abs(CumContentLeft-PercOfContent/2));

PercOfContent
PercOfLength=AxLeft(LeftIdx)+AxRight(RightIdx)


figure;
plot(AxLeft,CumContentLeft); hold on;
plot(AxRight,CumContentRight,'r-')