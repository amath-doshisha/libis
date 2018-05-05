function y=multi_cast(x,c)
if ~(nargin==2 && isa(c,'char') && (isa(x,'multi') || isfloat(x)))
    error('multi_cast(x,type)');
end
switch get_type(x)
    case 'r'
        switch c
            case 'r'
                y=multi('get_multi',x.data);
            case 'R'
                y=multi('get_imulti',x.data);
            case 'c'
                y=multi('complex',x.data);
            case 'C'
                y=multi('get_imulti',multi('complex',x.data).data);
            case 'd'
                y=double(x);
            case 'z'
                y=complex(double(x));
            otherwise
                error('multi_cast not supported: r->?');
        end
    case 'R'
        switch c
            case 'r'
                y=multi('get_multi',x.data);
            case 'R'
                y=multi('get_imulti',x.data);
            case 'c'
                y=multi('complex',multi('get_multi',x.data).data);
            case 'C'
                y=multi('complex',x.data);
            case 'd'
                y=double(x);
            case 'z'
                y=complex(double(x));
            otherwise
                error('multi_cast not supported: R->?');
        end
    case 'c'
        switch c
            case 'r'
                y=multi('real',x.data);
            case 'R'
                y=multi('real',multi('get_imulti',x.data).data);
            case 'c'
                y=multi('get_multi',x.data);
            case 'C'
                y=multi('get_imulti',x.data);
            case 'd'
                y=real(double(x));
            case 'z'
                y=double(x);
            otherwise
                error('multi_cast not supported: c->?');
        end
    case 'C'
        switch c
            case 'r'
                y=multi('real',multi('get_multi',x.data).data);
            case 'R'
                y=multi('real',x.data);
            case 'c'
                y=multi('get_multi',x.data);
            case 'C'
                y=multi('get_imulti',x.data);
            case 'd'
                y=real(double(x));
            case 'z'
                y=double(x);
            otherwise
                error('multi_cast not supported: C->?');
        end
    case 'd'
        switch c
            case 'r'
                y=multi(x);
            case 'R'
                y=imulti(x);
            case 'c'
                y=multi(complex(x));
            case 'C'
                y=imulti(complex(x));
            case 'd'
                y=x;
            case 'z'
                y=complex(x);
            otherwise
                error('multi_cast not supported: d->?');
        end
    case 'z'
        switch c
            case 'r'
                y=multi(real(x));
            case 'R'
                y=imulti(real(x));
            case 'c'
                y=multi(x);
            case 'C'
                y=imulti(x);
            case 'd'
                y=real(x);
            case 'z'
                y=x;
            otherwise
                error('multi_cast not supported: z->?');
        end
    otherwise
        error('multi_cast not supported: ?->?');
end
