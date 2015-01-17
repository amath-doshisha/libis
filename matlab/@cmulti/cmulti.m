classdef cmulti
    properties(SetAccess=private,GetAccess=private)
        data;
    end
    methods
        function obj=cmulti(cmd,prec,arg1,arg2)
            if nargin<=0
                obj.data=cmulti_allocate('default');
            elseif nargin==1
                obj.data=cmulti_allocate(cmd);   
            elseif nargin==2
                obj.data=cmulti_allocate(cmd,prec);   
            elseif nargin==3
                obj.data=cmulti_allocate(cmd,prec,arg1);
            else
                obj.data=cmulti_allocate(cmd,prec,arg1,arg2);                
            end
        end
    end
end
