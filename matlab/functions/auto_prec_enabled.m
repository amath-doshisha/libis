%% Enable auto prec mode.
function auto_prec_enabled()
global auto_prec_mode;
if exist('auto_prec_mode','var')==0 || length(auto_prec_mode)~=1
    auto_prec_mode=false;
end
auto_prec_mode=true;