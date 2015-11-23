% C=horzcat(A1,A2,...)
function C=horzcat(varargin)
cmd='horzcat';
i=1;
while i<=nargin
    if isa(varargin{i},'multi')
        varargin2{i}=varargin{i}.data;
    else
        varargin2{i}=multi(varargin{i}).data;
    end
    i=i+1;
end
C=multi(cmd,varargin2{:});