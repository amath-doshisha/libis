classdef multi
    properties(SetAccess=private,GetAccess=public)
        data;
    end
    methods
        function obj=multi(cmd,varargin)
            if nargin==1 && isa(cmd,'multi')
                obj.data=multi_mex('copy',get_default_prec(),double(get_auto_prec_mode()),cmd.data);
            elseif nargin==1 && isa(cmd,'double')
                obj.data=multi_mex('set_d',get_default_prec(),double(get_auto_prec_mode()),cmd);
            elseif nargin==1 && isa(cmd,'logical')
                obj.data=multi_mex('set_d',get_default_prec(),double(get_auto_prec_mode()),double(cmd));
            elseif nargin==1 && isa(cmd,'int32')
                obj.data=multi_mex('set_d',get_default_prec(),double(get_auto_prec_mode()),double(cmd));
            elseif nargin==1 && isa(cmd,'int64')
                obj.data=multi_mex('set_d',get_default_prec(),double(get_auto_prec_mode()),double(cmd));
            elseif nargin==1 && isa(cmd,'cell')
                obj.data=multi_mex('set_s',get_default_prec(),double(get_auto_prec_mode()),cmd);
            else
                obj.data=multi_mex(cmd,get_default_prec(),double(get_auto_prec_mode()),varargin{1:end});
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
%            A=multi('set_rand','r',varargin{1:end});
            A=multi(rand(varargin{1:end},'double'));
        end
        
        function A=eye(varargin)
            A=multi('set_eye','r',varargin{1:end});
        end
    end
end