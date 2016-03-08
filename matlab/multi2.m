function [out1,out2]=multi2(cmd,varargin)
out1=multi();
out2=multi();
if nargout<=1
    out1.data=multi_mex(cmd,get_default_prec(),double(get_auto_prec_mode()),varargin{1:end});
else
    [out1.data,out2.data]=multi_mex(cmd,get_default_prec(),double(get_auto_prec_mode()),varargin{1:end});
end
