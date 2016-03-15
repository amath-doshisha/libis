function [out1,out2]=multi2(cmd,varargin)
if nargout<=1
    out1=multi_mex(cmd,get_default_prec(),double(get_auto_prec_mode()),varargin{1:end});
else
    [out1,out2]=multi_mex(cmd,get_default_prec(),double(get_auto_prec_mode()),varargin{1:end});
end
