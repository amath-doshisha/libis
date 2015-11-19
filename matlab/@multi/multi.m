classdef multi
    properties(SetAccess=private,GetAccess=public)
        data;
    end
    methods
        function obj=multi(cmd,varargin)
            if nargin==1 && isa(cmd,'double')
                obj.data=multi_mex('set_d',get_default_prec(),cmd);
            elseif nargin==1 && isa(cmd,'multi')
                obj.data=multi_mex('copy',get_default_prec(),cmd.data);
            else
                obj.data=multi_mex(cmd,get_default_prec(),varargin{1:end});
            end
        end
    end
    
    methods(Static=true)
        function A=zeros(varargin)
            A=multi('set_zeros',varargin{1:end});
        end
        function A=ones(varargin)
            A=multi('set_ones',varargin{1:end});
        end
    end
end