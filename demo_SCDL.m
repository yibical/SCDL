close all;clear all;clc;
addpath('utilities');

addpath('.\data');

load AR_DR_DAT
%load UCF_SPORTED_Task1
%load YaleB_DR_DAT
tr_dat = Train_DAT;
tt_dat = Test_DAT;
trls = trainlabels;
ttls = testlabels;

%%%%%%%%%%%%%%%%%%%%%%%%
%SCDL parameter
%%%%%%%%%%%%%%%%%%%%%%%%
nClass        =  100;
nDCIter       =  15;
%show          =   true;
show          =   false;
lambda1 = 0.2; % regulization term
lambda2 = 0.2;    % parameter of ||Z_j^TZ_i||_F^2 in Eq.(1)
nDim = 300;  % feature dimension

%%%%%%%%%%%%%%%%%%%%%%%%
%SCDL input parameter
%%%%%%%%%%%%%%%%%%%%%%%%
opts.nClass        =   nClass;
opts.wayInit       =   'PCA';
opts.lambda1       =   lambda1;
opts.lambda2       =   lambda2;
opts.nDCIter         = nDCIter;
opts.show          =   show;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCDL initialize the dictionary and coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_ini = Eigenface_f(tr_dat,nDim);
tr_dat = P_ini'*tr_dat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCDL learn the dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[newDict,newDrls,newCoef,newCoeflabel]  = SCDL(tr_dat, trls, opts);
%%%%%%%%%%%%%%%%%%%%%%%%
% SCDL Classification
%%%%%%%%%%%%%%%%%%%%%%%%
tt_dat = P_ini'*tt_dat;
ID = SCDLSC(tt_dat,opts.nClass,newDict,newDrls,newCoef,newCoeflabel);

correct_rate = sum(ID'==ttls)/(length(ttls))
fprintf('%s%8f\n','reco_rate  =  ',correct_rate);
fid = fopen(['result\demo_SCDL_result_AR.txt'],'a');
fprintf(fid,'\n%s\n','==========================================');
fprintf(fid,'%s%8f\n','reco_rate1 = ',correct_rate);
fclose(fid);