function set_default_prec(prec)
global default_prec;
if exist('default_prec','var')==0 || length(default_prec)~=1
    default_prec=64;
end
if nargin<=0
    default_prec=64;
else
    if prec<=1
        warning('default_prec was set 2, due to the default_prec shoud be greater than or equal to 2.');
        default_prec=2;
    else
        default_prec=prec;
    end
end