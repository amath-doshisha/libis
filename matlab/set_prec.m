function set_prec(prec)
global default_prec;
if exist('default_prec','var')==0 || length(default_prec)~=1
    default_prec=64;
end
if nargin<=0
    default_prec=64;
else
    default_prec=prec;
end