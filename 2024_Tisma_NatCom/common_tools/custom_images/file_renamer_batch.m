function file_renamer_batch
%used to correct file names
workpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230129_BSG219-BSG217_comparison\';
namebase='BSG219_noXyl_003';
insert_text='xy1'
L_pre=length(namebase);
curpath=pwd; 
cd(workpath);
if ~isdir('tiffs_renamed'), mkdir('tiffs_renamed'); end
cd([workpath, 'tiffs\'])
filelist=dir('*.tif');
cd(curpath);
for ii=1:length(filelist);
    disp(length(filelist)-ii);
    this_tiffname=filelist(ii).name;
    prefix=this_tiffname(1:L_pre);
    suffix=this_tiffname(L_pre+1:end);
    thistiff=imread([workpath, 'tiffs\', this_tiffname]);
    new_tiffname=[prefix insert_text suffix];
    imwrite(thistiff,[workpath, 'tiffs_renamed\', new_tiffname]);
end
disp('done')