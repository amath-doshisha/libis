function set_default_round_mode(mode)
global default_round_mode;
if exist('default_round_mode','var')==0 || length(default_round_mode)~=1
    default_round_mode=0;
end
if nargin<=0
    default_round_mode=0;
else
    if mode>0
        default_round_mode=1;
    elseif mode<0
        default_round_mode=-1;
    else
        default_round_mode=0;
    end
end