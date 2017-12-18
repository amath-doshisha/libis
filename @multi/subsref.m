% y=subsref(x,s)
function varargout=subsref(x,s)
cmd='subsref';
switch s(1).type
    case '()'
        if ~isa(x,'multi')
            x=multi(x);
        end
        y=multi();
        y.data=subsref(x.data,s);
        varargout{1}=y;   
    otherwise
        varargout={builtin(cmd,x,s)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old version                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cmd='subsref';
% switch s(1).type
%     case '()'
%         if isa(x,'multi')
%             varargout{1}=multi(cmd,x.data,s);
%         else
%             varargout{1}=multi(cmd,multi(x).data,s);
%         end
%     otherwise
%         varargout={builtin(cmd,x,s)};
% end