function [Di_new,delet] = SCDL_UpdateDi (A,X,index,classn,SCDL_ipts,SCDL_par)

D    =  SCDL_ipts.D;
tau3 =  1;
tau2 =  1;
trls =  SCDL_ipts.trls;
dls =  SCDL_par.dls;
delet = [];

Do = D(:,dls~=index);
Di = D(:,dls==index);
Ai = A(:,trls==index);
Ao = A(:,trls~=index);
Xi = X(:,trls==index);
Xo = X(:,trls~=index);

Xi_i = Xi(dls==index,:);
Xi_o = Xi(dls~=index,:);
Xo_i = Xo(dls==index,:);
Xo_o = Xo(dls~=index,:);

Zi = Ai-Do*Xi_o;
Zo = Ao-Do*Xo_o;

for i = 1: size(Di,2)
    Yi = Zi - Di*Xi_i + Di(:,i)*Xi_i(i,:);
    Ui = Ai - Di*Xi_i + Di(:,i)*Xi_i(i,:);
    Vo = Zo - Di*Xo_i + Di(:,i)*Xo_i(i,:);
    
    UaXti = zeros(size(Yi*(Xi_i(i,:))'));
    for t_i = 1:classn
        if t_i~=index
        Xt_i = X(dls==index,trls==t_i);
        UaXti = UaXti + (0  - Di*Xt_i + Di(:,i)*Xt_i(i,:))*(Xt_i(i,:))';
        end
    end

    tem1 = -Yi*(Xi_i(i,:))' - (tau2)*Ui*(Xi_i(i,:))' - Vo*(Xo_i(i,:))' - tau3*UaXti;
    tem  = -tem1;
    if norm(tem,2)<1e-6
        Di(:,i) = zeros(size(tem));
        delet = [delet i];
    else
        Di(:,i) = tem./norm(tem,2);   
    end
end

Di_new = Di;