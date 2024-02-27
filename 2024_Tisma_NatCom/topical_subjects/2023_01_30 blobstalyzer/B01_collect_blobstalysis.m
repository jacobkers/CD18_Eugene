function B01_collect_blobstalysis
close all;
users=[{'Tisma'} {'Jacob'}];
collect_type='1st submission';
collect_type='rebuttal';
collect_type='final';

switch collect_type
    case '1st submission'  ,collect_first_submission(users);       
    case 'rebuttal', collect_rebuttal(users);
    case 'final', collect_final(users);    
end
