function [out1,out2,out3,out4]=eig_verify(A,varargin)
cmd='eig_verify';
if isa(A,'multi')
    if nargout<=1
        out1=multi();        
        out1.data=multi4(cmd,A.data,varargin{:});
    elseif nargout<=2
        out1=multi();        
        out2=multi();        
        [out1.data,out2.data]=multi4(cmd,A.data,varargin{:});
    elseif nargout<=3
        out1=multi();
        out2=multi();   
        out3=multi();   
        [out1.data,out2.data,out3.data]=multi4(cmd,A.data,varargin{:});
    elseif nargout<=4
        out1=multi();
        out2=multi();
        out3=multi();
        out4=multi();
        [out1.data,out2.data,out3.data,out4.data]=multi4(cmd,A.data,varargin{:});
    else
        error('multi/eig.m');
    end
else
    error('multi/eig.m');
end
