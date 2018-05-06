% y=subsref(x,s) i.e. y=x(s)
function varargout=subsref(x,s)
cmd='subsref';
switch s(1).type
    case '()'
        if isa(x,'multi')
            y=multi('empty',get_type(x));
            y.data=subsref(x.data,s);
            varargout{1}=y;
        else
            varargout={builtin(cmd,x,s)};
        end
    otherwise
        varargout={builtin(cmd,x,s)};
end