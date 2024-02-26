function [namestructure,L_pre]=neutralize_pre_part_of_name(thisname,pathinfo);
    %this function makes sure that letters in the name base will not
    %disturb the analysis of time , xy, c and z info    
    L_name=length(thisname);
    L_post=length([pathinfo.txycz_base, '.tif']);
    L_pre=L_name-L_post;
    real_post=thisname(L_pre+1:end);
    fake_pre=repmat('q', 1, L_pre);
    namestructure=[fake_pre, real_post];