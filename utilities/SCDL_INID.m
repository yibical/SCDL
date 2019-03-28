function D    =    SCDL_INID(data,nCol,wayInit)

m   =    size(data,1);

switch lower(wayInit)
    case {'ori'}
        D = data;
    case {'pca'}
        [D,disc_value,Mean_Image]   =    Eigenface_f(data,nCol-1);
        %[D,disc_value,Mean_Image]   =    Eigenface_f(data,7-1);
        D                           =    [D Mean_Image./norm(Mean_Image)];
    case {'random'}
        phi                         =    randn(m, nCol);
        phinorm                     =    sqrt(sum(phi.*phi, 2));
        D                           =    phi ./ repmat(phinorm, 1, nCol);
    otherwise 
        error{'Unkonw method.'}
end
return;
