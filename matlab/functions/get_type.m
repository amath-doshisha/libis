function c=get_type(x,y)
% 入力1ならそのタイプ 入力2なら2つを内包するタイプ
if nargin==1
    if ~isa(x,'multi')
        x=multi(x);
        if is_type_r(x)
            c='d';
        elseif is_type_c(x)
            c='z';
        else
            error('get_type: ???');
        end
    elseif isa(x,'multi')
        if  is_type_r(x)
            c='r';
        elseif is_type_c(x)
            c='c';
        elseif is_type_R(x)
            c='R';
        elseif is_type_C(x)
            c='C';
        else
            error('get_type: ???');
        end
    else
        error('get_type: ???');
    end
else
    cx=get_type(x);
    cy=get_type(y);
    if cx=='C' || cy=='C'
        c='C';
    elseif (cx=='R' && (cy=='d' || cy=='r' || cy=='R')) || (cy=='R' && (cx=='d' || cx=='r' || cx=='R'))
        c='R';
    elseif (cx=='R' && (cy=='z' || cy=='c')) || (cy=='R' && (cx=='z' || cx=='c'))
        c='C';
    elseif cx=='c' || cy=='c'
        c='c';
    elseif (cx=='r' && (cy=='d' || cy=='r' )) || (cy=='r' && (cx=='d' || cx=='r'))
        c='r';
    elseif (cx=='r' && cy=='z') || (cy=='r' && cx=='z')
        c='c';
    elseif cx=='z' || cy=='z'
        c='z';
    elseif (cx=='d' && cy=='d')
        c='d';
    else
        error('get_type: (?,?)');
    end
end