%% y=subsref(x,s)
function y=subsref(x,s)
cmd='subsref';
switch s(1).type
    case '()'
        if isa(x,'multi')
            y=multi(cmd,x.data,s);
        else
            y=builtin(cmd,x,s);
        end
    case '.'
        y=builtin(cmd,x,s);
    otherwise
        error('subsref_error');
end
