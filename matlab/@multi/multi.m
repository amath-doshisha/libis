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
    
    % initialize functions
    methods(Static=true)
        function A=zeros(varargin)
            A=multi('set_zeros','r',varargin{1:end});
        end
        
        function A=ones(varargin)
            A=multi('set_ones','r',varargin{1:end});
        end
        
        function A=rand(varargin)
            A=multi('set_rand','r',varargin{1:end});
        end
        
        function A=eye(varargin)
            A=multi('set_eye','r',varargin{1:end});
        end
    end
end