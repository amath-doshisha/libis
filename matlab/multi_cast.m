function y=multi_cast(x,c)
if ~(nargin==2 && isa(c,'char') && (isa(x,'multi') || isfloat(x)))
    error('multi_cast(x,type)');
end
cmd='get';
if isa(x,'multi')
    y=multi(cmd,x.data,c);
else
    y=multi(cmd,multi(x).data,c);
end
if ~isstruct(y.data)
    y=y.data;
end