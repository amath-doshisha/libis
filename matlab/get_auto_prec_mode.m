%% Get auto prec mode.
function mode=get_auto_prec_mode()
global auto_prec_mode;
if exist('auto_prec_mode','var')==0 || length(auto_prec_mode)~=1
    auto_prec_mode=false;
end
mode=auto_prec_mode;