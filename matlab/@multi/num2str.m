%% s=num2str(A,opt)
function s=num2str(A,opt)
cmd='num2str';
if isa(A,'multi')
    if nargin==1
        s=strjoincell(get_s(A));
    else
        s=strjoincell(get_s(A,opt));
    end
else
    if nargin==1
        s=num2str(A);
    else
        s=num2str(A,opt);
    end
end