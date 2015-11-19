classdef rmulti
    properties(SetAccess=private,GetAccess=public)
        data;
    end
    methods
        function obj=rmulti(arg1,arg2,arg3,arg4,arg5)
            prec=get_prec;
            if nargin==2 && arg1==-1
                obj.data=arg2;
            elseif nargin==1 && isa(arg1,'double')
                obj.data=rmulti_data(prec,'init','set',arg1);
            elseif nargin==4 && arg1==0
                obj.data=rmulti_data(prec,'init',arg2,arg3,arg4);
            elseif nargin==3 && arg1==1
                if isreal(arg3)
                    x=rmulti(arg3);
                else
                    x=arg3;
                end
                if isa(x,'rmulti')
                    obj.data=rmulti_data(prec,'one',arg2,x.data);
                else
                    error('rmulti, one');
                end
            elseif nargin==4 && arg1==2
                if isreal(arg3)
                    x=rmulti(arg3);
                else
                    x=arg3;
                end
                if isreal(arg4)
                    y=rmulti(arg4);
                else
                    y=arg4;
                end
                if isa(x,'rmulti') && isa(y,'rmulti')
                    obj.data=rmulti_data(prec,'two',arg2,x.data,y.data);
                else
                    error('rmulti, two');
                end
            else
                error('rmulti');
            end
        end
    end
    
    methods(Static=true)
        function A=zeros(m,n)
            A=rmulti(0,'zeros',m,n);
        end
        
        function A=ones(m,n)
            A=rmulti(0,'ones',m,n);
        end
        
        function A=eye(m,n)
            A=rmulti(0,'eye',m,n);
        end
        
        function A=rand(m,n)
            A=rmulti(0,'rand',m,n);
        end
    end
end