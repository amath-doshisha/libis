% x=subsasgn(x,s,y), i.e. x(s)=y
function x=subsasgn(x,s,varargin)
cmd='subsasgn';
switch s(1).type
    case '()'      
        if ~isa(x,'multi')
            x=multi(y);
        end
        y=varargin{1};
        if ~isa(y,'multi')
            y=multi(y);
        end
        x.data=subsasgn(x.data,s,y.data);
    otherwise
        x=builtin(cmd,x,s,varargin{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old version                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cmd='subsasgn';
% switch s(1).type
%     case '()'
%         y=varargin{1};
%         if isa(x,'multi') && isa(y,'multi')
%             x=multi(cmd,x.data,s,y.data);
%         elseif isa(x,'multi')
%             x=multi(cmd,x.data,s,multi(y).data);
%         elseif isa(varargin,'multi')
%             x=multi(cmd,multi(x).data,s,y.data);
%         else
%             x=multi(cmd,multi(x).data,s,multi(y).data);
%         end
%     otherwise
%         x=builtin(cmd,x,s,varargin{:});
% end