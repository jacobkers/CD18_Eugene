function data_smz=medsmooth(data,window,how)
halfspan=ceil(window/2);
le=length(data);
data_smz=zeros(le,1);
hs=0;
for i=1:le
    if i>halfspan & le-i>halfspan, hs=halfspan; end
    if i<halfspan  , hs=min([i-1,le-i]);end
    if le-i<halfspan, hs=min([i-1,le-i]); end
    switch how
        case 'Median', data_smz(i)=nanmedian(data(i-hs:i+hs));
        case 'Average',  data_smz(i)=nanmean(data(i-hs:i+hs));
    end
end