%% x=subsasgn(x,s,y), i.e. x(s)=y
function x=subsasgn(x,s,y)
cmd='subsasgn';
switch s(1).type
    case '()'
        if isa(x,'multi')
            if isa(y,'multi')
                x=multi(cmd,x.data,s,y.data);
            else
                x=multi(cmd,x.data,s,multi(y).data);
            end
        else
            x=builtin(cmd,x,s,y);
        end
    case '.'
        x=builtin(cmd,x,s,y);
    otherwise
        error('subsasgn_error');
end