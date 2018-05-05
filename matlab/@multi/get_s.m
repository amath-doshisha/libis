% s=get_s(A,opt)
function s=get_s(A,f,digits)
cmd='get_s';
if nargin<2
    f='g';
end
if nargin<3
    digits=5;
end
if ~isa(f,'char') || ~(strcmp(f,'f') || strcmp(f,'e') || strcmp(f,'g'))
    f='g';
end
if isa(A,'multi')
    s=multi(cmd,A.data,f,digits).data;
else
    s=multi(cmd,A,f,digits).data;
end