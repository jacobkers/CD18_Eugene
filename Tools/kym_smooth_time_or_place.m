function kymo=kym_smooth_time_or_place(kymo,span,how);
[tt,xx]=size(kymo);
switch how
    case 'time'
        for ii=1:xx
            kymo(:,ii)=smooth(kymo(:,ii),span,'moving');
        end
    case 'place'
        for ii=1:tt
            kymo(ii,:)=smooth(kymo(ii,:),span,'moving');
        end
end


