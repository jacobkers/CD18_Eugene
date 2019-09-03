function a_00000_code_info_autogenerator_Eugene
%-------------------------------------------------------------------
%Title(replace)
%
%Summary: Auto-generate file comments
%
%Approach
%This tool is intended to build an auto-updating quicksheet for use of custom code packages 
    % The program searches all m.files in a project directory for any
    %piece of text between two taglines (without the 'minus signs) :
    % 'Use-this-section-for-a-Quicksheet' 
    % 'End-of-Quicksheet-section'
    %These lines are used to indicate important sections to adjust settings
    %etc. or to find back places where data is saved.

    % Disclaimer: the associated custom code was developed for scientific use. 
    % For efficient use of developing time, limited time is spent on user interfaces. 
    % Therefore, the user is supposed to have reasonable knowledge of Matlab code 
    % and an understanding of the functions and algorithms, especially those
    % where the user is asked to adjust settings as below. Jacob Kerssemakers 2014 
%------------------------------------------------------------[JK14]
%
%Input
%
%Output
%
%References
%
%-------------------------------------------------------------------

%paste directories and your search string (no *wildcards*)
mainprojectdir='C:\Users\jkerssemakers\Dropbox\CD_recent\BN_CD18_Eugene\Matlab\';
projectdirlist={'2018_01_09 FlashTrack',...
                '2019_01_22 BoxTrack',...
                '2019_07_16 Condensin with Plectonemes',...
                '2019_31_12 Various Analyses',...
                'Matlab_tools_CD18EK'};
           
for pr=2:2 %length(projectdirlist)
    clc
    disp('AUTO-GENERATED CODE DESCRIPTION');
    disp('Jacob Kerssemakers')
    disp(date);
    projectname=projectdirlist{pr};
    matlabsubpaths={''};   
    singledir=[mainprojectdir, projectname, '\']; 
    do_singledir=0;  %Use this if you want to search one directory, no subdirectories

    dt = datestr(now, 'mmmm_dd_yyyy_HH_MM');
    nm=strcat(mainprojectdir,'a0_code_info_',projectname,'_version_',dt); 
    diary(strcat(nm,'.txt'));    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('__________________________________________________________');
    disp('Order of listing:');
    disp('Main programs: project-specific shell programs');  %A
    disp('Secondary programs are project-specific sub-programs'); %B 
    disp('"Tools are small sub-functions'); %C
    disp('"Favorite Tools" are handy sub-functions of common interest'); %C*
    disp('_________________________________________________________');
    disp('       ');
    Headers=[{'A'} , {'B'} {'C'} {'C*'}];  
    display_options=[{'Function list:'} , {'Full description:'}];  
    %display_options=[{'Function list:'}];  
    %Headers=[{'A'},{'B'}];  
    %these terms label comments; to be searched in that order
    NHeaders=length(Headers);

    progdir=pwd;
    Nsubdirs=length(matlabsubpaths);

    if do_singledir, lp=1;end
    for dd=1:length(display_options)
        display_option=display_options{dd};
        disp(' ');
        disp(display_option);
        totalcodelines=0;
        totalcharcount=0;
        totalfunctioncount=0;
        foundmatches=0;
        for hh=1:NHeaders
            HeaderLabel=char(Headers(hh));
            switch HeaderLabel
                case 'A', titleheader='Code category: main analysis';    
                case 'B', titleheader='Code category: sub analysis';    
                case 'C', titleheader='Code category: tools' ; 
                case 'C*', titleheader='Code category: favorite tools' ;   
            end    
            disp(' ');
            disp(titleheader);
            disp('________________');

            for ii=1:Nsubdirs %for all subdirpaths
                p = genpath(char(strcat(singledir,matlabsubpaths(ii)))); 
                %String containing all (sub)dirs
                lp=length(p);
                cnt=0;
                oldi=1; nwi=oldi;
                while cnt<lp                        %separate out directory names
                    cnt=cnt+1;
                    nwi=nwi+1;
                    str=char(p(cnt));
                    if strcmp(str,';')              %This is one directory name
                        dirnm=p(oldi:nwi-2);
                        if do_singledir, 
                        dirnm=singledir;
                        cnt=lp;
                        end
                        dirnm;
                        cd(dirnm);                       %read the m.file names in this directory
                        FilNames=dir('*.m');
                        lm=length(FilNames);
                        for j=1:lm  %for all matlab files found
                            nm=FilNames(j).name;
                            xx=strfind(nm, '._');
                            if length(xx)<1
                                filecontent=textread(nm, '%s', 'delimiter', '\n','whitespace', '');                   
                                mfilestruct=Evaluate_m_File(filecontent,HeaderLabel);
                                totalcodelines=totalcodelines+mfilestruct.mfilelength;
                                totalcharcount=totalcharcount+mfilestruct.charactercount;
                                totalfunctioncount=totalfunctioncount+mfilestruct.functioncount;
                                foundmatches=foundmatches+mfilestruct.tagline_number;
                                Evaluate_functions_in_mfile(filecontent,mfilestruct,titleheader,display_option);
                            end
                        end
                        oldi=nwi;
                        cd(progdir);
                    end     
                end
            end
        end
    end

    
    disp(' ');
    disp(strcat('Matlab code associated with project:',singledir));

    %disp(' ');
    disp(strcat('number of characters:',num2str(totalcharcount)));
    disp(strcat('number of code lines:',num2str(totalcodelines)));
    disp(strcat('number of functions:',num2str(totalfunctioncount)));
    disp(strcat('number of comment matches:',num2str(foundmatches)));
    disp('__________________________________________________________');
    diary off
end

function mfilestruct=Evaluate_m_File(content,HeaderLabel)
%This function evaluates a piece of matlab code 'content' for
%matlab functions, comments etc. It is intended to allow the user to see if
%functions are of interest to store separately.
   
    mfilestruct.tagline_number=0;        %number of taglines found
    mfilestruct.tagline_start=[];        %sets start line nrs. of sections enclosed by taglines
    mfilestruct.tagline_end=[];          %sets end line nrs. of sections enclosed by taglines
    lc=length(content);
    mfilestruct.mfilelength=lc;
    mfilestruct.charactercount=0;
    mfilestruct.functioncount=0;
    for k=1:lc  %read in all lines
        lne=(char(content(k)));  %this is the current line 
        mfilestruct.charactercount=mfilestruct.charactercount+length(lne);
        %various text to recognize comments
       StartLabel=strcat('JWJK_',HeaderLabel,':');
       StopLabel=strcat(':','JWJK_',HeaderLabel);
       instructstart=strfind(lne, StartLabel);              
       instructstop=strfind(lne, StopLabel);        
       
        foundfunction=strfind(lne,'function');
        commentline=strfind(lne,'%');
        if length(foundfunction)>0 & isempty(commentline)
            mfilestruct.functioncount=mfilestruct.functioncount+1;
        end
        if length(instructstart)>0  %found the sentence   
                mfilestruct.tagline_number=mfilestruct.tagline_number+1;
                mfilestruct.tagline_start=[mfilestruct.tagline_start k]; 
        end  
        if length(instructstop)>0  %found the sentence             
            %mfilestruct.tagline_number=mfilestruct.taglinenumber+1;
            mfilestruct.tagline_end=[mfilestruct.tagline_end k]; 
        end  
    end
    mfilestruct.tagline_start=[mfilestruct.tagline_start lc];
    mfilestruct.tagline_end=[mfilestruct.tagline_end lc];
   
    
    
function Evaluate_functions_in_mfile(filecontent,mfilestruct,titleheader,display_option);
%This function takes function-headered sections of code out of a file; the
%user can decide to watch it and to judge its uefullness. JacobKers2014
%mfilestruct.functionnumber       %number of functions found
%mfilestruct.functionheaderlines  %line indices, including the last line
flines=mfilestruct.tagline_start;
flines_end=mfilestruct.tagline_end;
f=flines(end);
for i=1:length(flines)-1;
    switch display_option
        case 'Full description:'            
            msg0=char(filecontent(flines(1)-1));
            [outparams, functionname, inparams]=split_function_call(msg0); 
            disp(['Name: ', functionname]);
            disp(['Parameters in: ' ,inparams]);
            disp(['Parameters out: ' ,outparams]);
            
            disp(titleheader);            
            for k=flines(i)+1:flines_end(i)-1  %Display the instruction lines
              msg3=char(filecontent(k));
              prc = strfind(msg3, '%');
              if ~isempty(k), msg3(prc)=[]; end
              disp(msg3);
            end
            [~,sourcename,~]=split_function_call(char(filecontent(1)))
            msg1=strcat('Source code: ',sourcename, ' line:',num2str(flines(i)),' '); 
            disp(msg1);
            disp('__________________________________________________________ ');         
            for jj=1:1,disp(' ');end
            case 'Function list:'
                msg0=strcat(char(filecontent(flines(1)-1)));  %function call line
                [outparams, functionname, inparams]=split_function_call(msg0);        
                disp(functionname);
    end       
    end  

 function [outparams, functionname, inparams]=split_function_call(calltext);
 %split the function call
 idx0=strfind(calltext, 'function'); %first text after word function
 idx1=strfind(calltext, '=');
 
 if ~isempty(idx1)
    outparams=calltext(idx0+8:idx1-1);
    idx2a=strfind(outparams, '[');
    idx2b=strfind(outparams, ']');
    if ~isempty(idx2a) %mutiple output
      outparams=outparams(idx2a+1:idx2b-1);     
    end  %no input
 else
     outparams='no output';
 end
 
 idx3a=strfind(calltext, '(');
 idx3b=strfind(calltext, ')');
 
 if ~isempty(idx3a)
      functionname=calltext(idx1+1:idx3a-1);
      inparams=calltext(idx3a+1:idx3b-1);
 else  %no input
     if ~isempty(idx1)  
              functionname=calltext(idx1+1:end);
     else %no output
         functionname=calltext(idx0+9:end);
     end
     inparams='no input';
 end
 
 dum=1;
    
function Manual_Template_Kers2017
%JWJK_:----[add ABCorC*----------------------------------------------------
%Title(replace)
%Summary
%Approach
%Input
%Output
%References
%:JWJK_-----[add ABCorC*---------------------------------------------------


