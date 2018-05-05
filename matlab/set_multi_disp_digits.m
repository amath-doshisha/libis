function set_multi_disp_digits(digits)
global multi_disp_digits;
if exist('multi_disp_digits','var')==0 || length(multi_disp_digits)~=1
    multi_disp_digits=15;
end
if nargin<=0
    multi_disp_digits=15;
else
    if digits<0
        multi_disp_digits=0;
    else
        multi_disp_digits=digits;
    end
end