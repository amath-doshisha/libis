function y=multi_cast(x,c)

if nargin==2
    if ~isa(c,'char')
        error('multi_cast(multi or (d or z),char) needed');
    end
    if isa(x,'multi') || isfloat(x) % multi ‚© single ‚© double
        switch get_type(x)
            case 'r'
                switch c
                    case 'r'
                        y=multi('copy',x.data);
                    case 'R'
                        y=multi('icopy',x.data);
                    case 'c'
                        y=multi('complex',x.data);
                    case 'C'
                        y=multi('icopy',multi('complex',x.data).data);
                    case 'd'
                        y=double(x);
                    case 'z'
                        y=double(x);
                    otherwise
                        error('multi_cast not supported:(r,?)');
                end
            case 'R'
                switch c
                    case 'r'
                        y=multi('copy',x.data);
                    case 'R'
                        y=multi('icopy',x.data);
                    case 'c'
                        y=multi('complex',multi('copy',x.data).data);
                    case 'C'
                        y=multi('complex',x.data);
                    case 'd'
                        y=double(x);
                    case 'z'
                        y=double(x);
                    otherwise
                        error('multi_cast not supported:(R,?)');
                end
            case 'c'
                switch c
                    case 'c'
                        y=multi('copy',x.data);
                    case 'C'
                        y=multi('icopy',x.data);
                    case 'z'
                        y=double(x);
                    otherwise
                        error('multi_cast not supported:(c,?)');
                end
            case 'C'
                switch c
                    case 'c'
                        y=multi('copy',x.data);
                    case 'C'
                        y=multi('icopy',x.data);
                    case 'z'
                        y=double(x);
                    otherwise
                        error('multi_cast not supported:(C,?)');
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
                        error('multi_cast not supported:(d,?)');
                end
            case 'z'
                switch c
                    case 'c'
                        y=multi(x);
                    case 'C'
                        y=imulti(x);
                    case 'z'
                        y=x;
                    otherwise
                        error('multi_cast not supported:(z,?)');
                end
            otherwise
                error('multi_cast not supported:(?,?)');
        end
    else
        error('multi_cast(multi or (d or z),char) needed');
    end
else
    error('multi_cast(multi or (d or z),char) needed');
end