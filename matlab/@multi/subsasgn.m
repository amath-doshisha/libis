% x=subsasgn(x,s,y), i.e. x(s)=y
function x=subsasgn(x,s,varargin)
cmd='subsasgn';
switch s(1).type
    case '()'
        y=varargin{1};
        if ~isa(x,'multi')
            x=multi(x);
        end
        if ~isa(y,'multi')
            y=multi(y);
        end
        cx=get_type(x);
        cy=get_type(y);
        c=get_type(x,y);
        if cx~=c
            x=multi_cast(x,c);
        end
        if cy~=c
            y=multi_cast(y,c);
        end
        size_old=size(x.data);
        x.data=subsasgn(x.data,s,y.data);
        size_new=size(x.data);
%        if ~isequal(size_old,size_new)
%            error('The size is changed.');
%        end
    otherwise
        x=builtin(cmd,x,s,varargin{:});
end