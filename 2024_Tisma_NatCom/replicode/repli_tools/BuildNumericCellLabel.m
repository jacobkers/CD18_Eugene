function NumCellLabel=BuildNumericCellLabel(cellno);
    %function transforms label to unique number. For movies, three last
    %decimals are reserved for frame numbers
        NumCellLabel=str2num(cellno);  %regular all-number label
        if isempty(NumCellLabel) %label conains string at end
            [numpart,~, ~, nextindex]=sscanf(cellno,'%d');
             numpart2=str2num(cellno(nextindex+2:end));
             NumCellLabel=1000*numpart+numpart2;
             if isempty(NumCellLabel)  %still not good...
                 sel=find(~isletter(cellno));
                 NumCellLabel=str2num([cellno(sel(2:end))])+100000;
             end
        end 
        dum=1;
    