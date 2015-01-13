classdef rmulti
    properties(SetAccess=private,GetAccess=private)
        data;
    end
    methods
        function obj=rmulti(cmd,prec,arg1,arg2)
            if nargin<=0
                obj.data=rmulti_allocate('default');
            elseif nargin==1
                obj.data=rmulti_allocate(cmd);   
            elseif nargin==2
                obj.data=rmulti_allocate(cmd,prec);   
            elseif nargin==3
                obj.data=rmulti_allocate(cmd,prec,arg1);
            else
                obj.data=rmulti_allocate(cmd,prec,arg1,arg2);                
            end
        end
    end
end