function [out1,out2,out3,out4,out5]=multi5(cmd,varargin)
if nargout<=1
    out1=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
elseif nargout<=2
    [out1,out2]=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
elseif nargout<=3
    [out1,out2,out3]=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
elseif nargout<=4
    [out1,out2,out3,out4]=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
elseif nargout<=5
    [out1,out2,out3,out4,out5]=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
else
    error('multi5: Too many outputs.');
end
