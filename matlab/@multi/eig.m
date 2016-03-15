function [out1,out2]=multi_eig(A,varargin)
cmd='eig';
if isa(A,'multi')
    if nargout<=1
        out1=multi();        
        out1.data=multi2(cmd,A.data,varargin{:});
    else
        out1=multi();        
        out2=multi();        
        [out1.data,out2.data]=multi2(cmd,A.data,varargin{:});
    end
else
    error('multi/eig.m');
end
