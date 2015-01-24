classdef cmulti
    properties(SetAccess=private,GetAccess=public)
        data;
    end
    methods
        function obj=cmulti(arg1,arg2,arg3,arg4,arg5)
            prec=get_prec;
            if nargin==2 && arg1==-1
                obj.data=arg2;
            elseif nargin==1 && isa(arg1,'double')
                obj.data=cmulti_data(prec,'init','set',arg1);
            elseif nargin==1 && isa(arg1,'rmulti')
                obj.data=cmulti_data(prec,'init','set_r',arg1.data);
            elseif nargin==4 && arg1==0
                obj.data=cmulti_data(prec,'init',arg2,arg3,arg4);
            elseif nargin==3 && arg1==1
                if isa(arg3,'double')
                    x=cmulti(arg3);
                elseif isa(arg3,'rmulti')
                    x=cmulti(arg3);
                else
                    x=arg3;
                end
                if isa(x,'cmulti')
                    obj.data=cmulti_data(prec,'one',arg2,x.data);
                else
                    error('cmulti, one');
                end
            elseif nargin==4 && arg1==2
                if isa(arg3,'double')
                    x=cmulti(arg3);
                elseif isa(arg3,'rmulti')
                    x=cmulti(arg3);
                else
                    x=arg3;
                end
                if isa(arg4,'double')
                    y=rmulti(arg4);
                elseif isa(arg4,'rmulti')
                    y=cmulti(arg4);
                else
                    y=arg4;
                end
                if isa(x,'cmulti') && isa(y,'cmulti')
                    obj.data=cmulti_data(prec,'two',arg2,x.data,y.data);
                else
                    error('cmulti, two');
                end
            else
                error('cmulti');
            end
        end
    end
    
    methods(Static=true)
        function A=zeros(m,n)
            A=cmulti(0,'zeros',m,n);
        end
        
        function A=ones(m,n)
            A=cmulti(0,'ones',m,n);
        end
        
        function A=eye(m,n)
            A=cmulti(0,'eye',m,n);
        end
        
        function A=rand(m,n)
            A=cmulti(0,'rand',m,n);
        end
    end
end
