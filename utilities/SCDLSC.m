function [ID] = SCDLSC(TestDat,NumClass,Dict,Drls,Coef, CoefL)

%%%%%%%%%%%%%%%%%%
% normalize energy
%%%%%%%%%%%%%%%%%%

TestDat = TestDat*diag(1./sqrt(sum(TestDat.*TestDat)));

% Add modified code 
Coef = Coef*diag(1./sqrt(sum(Coef.*Coef)));
% Modified end

lambda   =   0.005;

weight   =  5.6;

nClass   =   NumClass;

td1_ipts.D    =   Dict;
td1_ipts.tau1 =   lambda;
if size(td1_ipts.D,1)>=size(td1_ipts.D,2)
   td1_par.eigenv = eigs(td1_ipts.D'*td1_ipts.D,1);
else
   td1_par.eigenv = eigs(td1_ipts.D*td1_ipts.D',1);  
end

ID   =   [];
strength = [];
nNum     =  size(Dict,2);

Proj_M   =  inv(Dict'*Dict+lambda*eye(nNum))*Dict';
for indTest = 1:size(TestDat,2)
    %fprintf(['Totalnum:' num2str(size(TestDat,2)) 'Nowprocess:' num2str(indTest) '\n']);
    td1_ipts.y          =      TestDat(:,indTest);   
    %[opts]              =      IPM_SC(td1_ipts,td1_par);
    %s                   =      opts.x;
    s=Proj_M*(td1_ipts.y);

    

    
    for indClass  =  1:nClass
        temp_s            =  zeros(size(s));
        temp_s(indClass==Drls) = s(indClass==Drls);
        zz                =  TestDat(:,indTest)-td1_ipts.D*temp_s;
        gap(indClass)     =  zz(:)'*zz(:);
        
        temp_Coef = Coef;
        temp_Coef(:,CoefL == indClass) = 0;
        indClassNum = sum(CoefL == indClass);
        gCoef3(indClass)    =  (norm(s'*temp_Coef)^2)*(weight/indClassNum);
    end
    wgap3  = gap +gCoef3;
    index3 = find(wgap3==min(wgap3));
    id3    = index3(1);
    ID     = [ID id3];
    strength = [strength wgap3];
end  
ID = ID';
strength = strength';