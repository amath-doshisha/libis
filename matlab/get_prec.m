function prec=get_prec()
global default_prec;
if exist('default_prec','var')==0 || length(default_prec)~=1
    default_prec=64;
end
prec=default_prec;