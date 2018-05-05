classdef multi
    %    properties(SetAccess=private,GetAccess=public)
    %    properties
    properties(SetAccess=public,GetAccess=public)
        data;
    end
    methods
        function obj=multi(cmd,varargin)
            if nargin==0
                obj.data=multi_mex('set_zeros',get_default_prec(),double(get_default_round_mode()),'r',0,0);
            elseif nargin==1 && isa(cmd,'multi')
                obj.data=multi_mex('get_multi',get_default_prec(),double(get_default_round_mode()),cmd.data);
            elseif nargin==1 && isa(cmd,'double')
                obj.data=multi_mex('get_multi',get_default_prec(),double(get_default_round_mode()),cmd);
            elseif nargin==1 && (isnumeric(cmd) || islogical(cmd))
                obj.data=multi_mex('get_multi',get_default_prec(),double(get_default_round_mode()),double(cmd));
            elseif nargin==1 && isa(cmd,'cell')
                obj.data=multi_mex('get_multi',get_default_prec(),double(get_default_round_mode()),cmd);
            elseif nargin==1 && isa(cmd,'char')
                obj.data=multi_mex('get_multi',get_default_prec(),double(get_default_round_mode()),strsplitcell(cmd));
            else
                obj.data=multi_mex(cmd,get_default_prec(),double(get_default_round_mode()),varargin{1:end});
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
        
        function A=nan(varargin)
            A=multi('set_nan','r',varargin{1:end});
        end
        
        function A=inf(varargin)
            A=multi('set_inf','r',varargin{1:end});
        end
        
        
        function A=rand(varargin)
            %A=multi('set_rand','r',varargin{1:end});
            A=multi(rand(varargin{1:end},'double'));
        end
        
        function A=eye(varargin)
            A=multi('set_eye','r',varargin{1:end});
        end
    end
end