function A001_Erase_and_Build_ResultDirs(initval)
%JWJK_A:-------------------------------------------------------------------
%Description: For first time use of experiment; Check well before you use it!!

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------


if nargin<1
    usr='Jacob', batchrunindex=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end

%yesno = input(strcat('overwrite:', initval.pth_repli, '? y/n '),'s');
yesno='y';
if strcmp(yesno,'y')
    disp('overwriting...')
if 0 %isdir(initval.pth_repli), 
    rmdir(initval.pth_repli,'s'); 
end
    mkdir(initval.pth_repli);
    MatFilePath=strcat(initval.pth_repli,'ResultsPerCellMatlab\',initval.DirSep);
    mkdir(MatFilePath);
end

