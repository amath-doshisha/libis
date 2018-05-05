function digits=get_multi_disp_digits()
global multi_disp_digits;
if exist('multi_disp_digits','var')==0 || length(multi_disp_digits)~=1
    multi_disp_digits=15;
end
digits=multi_disp_digits;