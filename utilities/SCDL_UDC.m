function [Dict,Drls, Coef, Coeflabel] = SCDL_UDC(TrainDat, TrainLabel, opts, Dict_ini, Dlabel_ini, Coef_ini)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%
% normalize energy
%%%%%%%%%%%%%%%%%%
TrainDat = TrainDat*diag(1./sqrt(sum(TrainDat.*TrainDat)));

 NDL_par.dls        =     Dlabel_ini;
 NDL_ipts.D         =     Dict_ini;
 NDL_ipts.trls      =     TrainLabel;
 NDL_par.tau        =     opts.lambda1;
 
 
 lambda_a           =     opts.lambda1;
 NDL_par.lambda2    =     opts.lambda2;

 
 lambda_b           =     opts.lambda2;
 
 NDL_nit            =     1;
 drls                =    Dlabel_ini;
 coef                =     Coef_ini;
 while NDL_nit<=opts.nDCIter  
    if size(NDL_ipts.D,1)>size(NDL_ipts.D,2)
      NDL_par.c        =    1.05*eigs(NDL_ipts.D'*NDL_ipts.D,1);
    else
      NDL_par.c        =    1.05*eigs(NDL_ipts.D*NDL_ipts.D',1);
    end
    %-------------------------
    %updating the coefficient
    %-------------------------
    
    Dict_totsum = zeros(size(NDL_ipts.D,2));
    for ci = 1:opts.nClass
        Tmp_Dict_Part=NDL_ipts.D;
        Tmp_Dict_Part(:,NDL_par.dls ~= ci) = 0;
        Dict_totsum = Dict_totsum+Tmp_Dict_Part'*Tmp_Dict_Part;
    end
    
    for ci = 1:opts.nClass
        fprintf(['Updating coefficients, class: ' num2str(ci) '\n'])
        NDL_ipts.X         =  TrainDat(:,TrainLabel==ci);
        
        Tmp_Dict_Part=NDL_ipts.D;
        Tmp_Dict_Part(:,NDL_par.dls ~= ci) = 0;
        Tmp_Coef_Part = coef;
        Tmp_Coef_Part(:, TrainLabel ~= ci) = 0;
        Dict_sin = Tmp_Dict_Part;
        Coef_sin = Tmp_Coef_Part;
        PinvD = inv(NDL_ipts.D'*NDL_ipts.D+Dict_totsum+ lambda_b*(Coef_sin*Coef_sin')+lambda_a*eye(size(NDL_ipts.D,2)));
        afa = PinvD*((NDL_ipts.D'+Dict_sin')*NDL_ipts.X);
        
        Copts.A = afa;
        coef(:,TrainLabel==ci) =    Copts.A;    
    end
    
    %------------------------
    %updating the dictionary
    %------------------------
    for ci = 1:opts.nClass
        fprintf(['Updating dictionary, class: ' num2str(ci) '\n']);     
        [NDL_ipts.D(:,drls==ci),Delt(ci).delet]= SCDL_UpdateDi (TrainDat,coef,...
            ci,opts.nClass,NDL_ipts,NDL_par);
    end    
    newD = []; newdrls = []; newcoef = [];
    for ci = 1:opts.nClass
        delet = Delt(ci).delet;
        if isempty(delet)
           classD = NDL_ipts.D(:,drls==ci);
           newD = [newD classD];
           newdrls = [newdrls repmat(ci,[1 size(classD,2)])];
           newcoef = [newcoef; coef(drls==ci,:)];
        else
            temp = NDL_ipts.D(:,drls==ci);
            temp_coef = coef(drls==ci,:);
            for temp_i = 1:size(temp,2)
                if sum(delet==temp_i)==0
                    newD = [newD temp(:,temp_i)];
                    newdrls = [newdrls ci];
                    newcoef = [newcoef;temp_coef(temp_i,:)];
                end
            end
        end
    end
    
    NDL_ipts.D  = newD;
    coef         = newcoef;
    drls         = newdrls;
    NDL_par.dls        =     drls;
    
    NDL_nit = NDL_nit +1;
    end
    
    Dict = NDL_ipts.D;
    Drls = drls;

    % New modified code
    Coef = coef;
    Coeflabel = TrainLabel;
    % New modified end
return;