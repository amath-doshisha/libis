function y=imulti(x0,x1)
cmd='get_imulti';
if nargin==1
    if isa(x0,'double')
        y=multi(cmd,x0);
    elseif isnumeric(x0) || islogical(x0)
        y=multi(cmd,double(x0));
    elseif isa(x0,'multi')
        y=multi(cmd,x0.data);
    else
        y=multi(cmd,multi(x0).data);
    end
elseif nargin==2
    if isa(x0,'double') && isa(x1,'double') 
        y=multi(cmd,x0,x1);
    elseif (isnumeric(x0) || islogical(x0)) && (isnumeric(x1) || islogical(x1))
        y=multi(cmd,double(x0),double(x1));
    elseif (isa(x0,'multi') && isa(x1,'multi'))
        y=multi(cmd,x0.data,x1.data);
    else
        error('imulti: not supported');
    end
else
    error('imulti: nargin');
end