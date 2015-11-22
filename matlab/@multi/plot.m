function plot(varargin)
i=1;
while i<=nargin
    if isa(varargin{i},'multi')
        varargin2{i}=double(varargin{i});
    else
        varargin2{i}=varargin{i};
    end
    i=i+1;
end
plot(varargin2{:});