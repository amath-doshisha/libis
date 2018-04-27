function mode=get_default_round_mode()
global default_round_mode;
if exist('default_round_mode','var')==0 || length(default_round_mode)~=1
    default_round_mode=0;
end
mode=default_round_mode;