function [ Dict,Drls, Coef, Coeflabel] = SCDL(tr_dat, trls, opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%
% normalize energy
%%%%%%%%%%%%%%%%%%
Tr_Dat = tr_dat*diag(1./sqrt(sum(tr_dat.*tr_dat)));

%%%%%%%%%%%%%%%%%%
%initialize dict
%%%%%%%%%%%%%%%%%%
Dict_ini  =  []; 
Dlabel_ini = [];
for ci = 1:opts.nClass
    cdat          =    Tr_Dat(:,trls==ci);
    dict          =    SCDL_INID(cdat,size(cdat,2),opts.wayInit);
    Dict_ini      =    [Dict_ini dict];
    Dlabel_ini    =    [Dlabel_ini repmat(ci,[1 size(dict,2)])];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize coef of SCDL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ini_par.tau         =     opts.lambda1;
lambda_a            =     opts.lambda1;
ini_par.lambda      =     opts.lambda2;
lambda_b            =     opts.lambda2;

ini_ipts.D          =     Dict_ini;
ini_ipts.Dlabel     =     Dlabel_ini;

coef = zeros(size(Dict_ini,2),size(Tr_Dat,2));
if size(Dict_ini,1)>size(Dict_ini,2)
      ini_par.c        =    1.05*eigs(Dict_ini'*Dict_ini,1);
else
      ini_par.c        =    1.05*eigs(Dict_ini*Dict_ini',1);
end

Dict_totsum = zeros(size(ini_ipts.D,2));
for ci = 1:opts.nClass
    Tmp_Dict_Part=ini_ipts.D;
    Tmp_Dict_Part(:,ini_ipts.Dlabel ~= ci) = 0;
    Dict_totsum = Dict_totsum+Tmp_Dict_Part'*Tmp_Dict_Part;
end
for ci =  1:opts.nClass
    fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
    ini_ipts.X      =    Tr_Dat(:,trls==ci);
    Tmp_Dict_Part=ini_ipts.D;
    Tmp_Dict_Part(:,ini_ipts.Dlabel ~= ci) = 0;
    Dict_sin = Tmp_Dict_Part;
    ini_opts.A = inv(ini_ipts.D'*ini_ipts.D+Dict_totsum+(lambda_a+lambda_b)*eye(size(ini_ipts.D,2)))*((ini_ipts.D'+Dict_sin')*ini_ipts.X);
    coef(:,trls ==ci) =    ini_opts.A;
end

Coef_ini = coef;

[newDict,newDrls, newCoef, newCoeflabel]=SCDL_UDC(tr_dat, trls, opts, Dict_ini, Dlabel_ini, Coef_ini);
Dict= newDict;
Drls = newDrls;
Coef = newCoef;
Coeflabel=  newCoeflabel;

return;

