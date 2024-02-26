function [ind,digits,strtext,strcode,pre_str]=Get_TimePlaceIndicesFromName(filnamefull,str);
    filname=filnamefull(1:end-4);
    Lstr=length(str);
    numbrs=~isletter(filname);
    ind = max(strfind(filname,str));  %last 'str'
    if ~isempty(ind) 
     %if the string is found, count the numbers behind it 
        digits=0; keepcounting=1; counter=0; idx=0;
        while (keepcounting&idx<length(numbrs))
            idx=ind+Lstr+counter;
            if numbrs(idx)&(idx<=length(numbrs))
                digits=digits+1;
            else
                keepcounting=0;
            end
            counter=counter+1;
        end
        if digits>0  
            strtext=filname(ind:ind+Lstr+digits-1);
            strcode=strcat('%0',num2str(digits),'i');
        else      %false alarm     
            ind=[]; digits=[];  strtext=[];   strcode=[];
        end
    else
            ind=[]; digits=[];  strtext=[];   strcode=[];
    end
    if ind>1
        pre_str=filnamefull(1:ind-1);
    else
        pre_str=[];
    end
   