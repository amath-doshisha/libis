clc
%format long

AR=[1 1 3; 2 0 4; -1 6 -1];
BR=[2 19 8];

AC=[1+i 1-i 3+2*i; 2 i 4; -1 6 -1+i];
BC=[2+i 19-i 8+i];

ar=multi(AR);
br=multi(BR);

ac=multi(AC);
bc=multi(BC);

r=1;
c=1+i;

%{
%丸めの確認
format long
double(double(br\ar))
BR\AR
double(BR\AR)
%}

disp('start');
%{
%丸めが発生しているため、不可。
if isequal(double(br\ar),BR\AR)==false disp('false');end
if isequal(double(bc\ar),BC\AR)==false disp('false');end
if isequal(double(br\ac),BR\AC)==false disp('false');end
if isequal(double(bc\ac),BC\AC)==false disp('false');end
%}

%{
%この演算では、一番大きい値について、除算を行い、それ以外の要素については、零埋めする。
[2 19 8]\1

[2 2 19 8]\1
[2 19 8 8]\1

[2 2 19 8 8]\1
[2 2 19 8 8 8]\1

[2 2 256 8 8 8]\1
[256 2 256.1 8 8 8]\1

[1 1 5 1 1 5]\1

[1 -1 1 1 1 1]\1
BR\r%例外処理!!!
%}

%{
%この演算では、一番大きい値について、除算を行い、それ以外の要素については、零埋めする。

%実部の絶対値を比較、等しい場合は、さらに虚部の絶対値を比較し、大きい方の値をとる。
%全く同じ値の場合は、先に出現した値を見る。

[1+i 1 1 1 1 1]\1
[1-i 1 1 1 1 1]\1
[1-i 1+i 1 1 1 1]\1
[1+i 1-i 1 1 1 1]\1

[1+2i 1+4i 1 1 1 1]\1
[1+2i 1-4i 1 1 1 1]\1
%}
%{
[10+10i 1-100i 1 1 1 1]\1
1/(10+10i)
abs(1/(10+10i))
1/(10-100i)
abs(1/(10-100i))

[10+10i 10+11i 10-11i 9+11i 1 1]\1
1/(10+10i)
abs(1/(10+10i))
1/(10+11i)
abs(1/(10+11i))
1/(10-11i)
abs(1/(10-11i))
%}


if isequal(double(br\r),BR\r)==false disp('false');end
if isequal(double(br\c),BR\c)==false disp('false');end
if isequal(double(bc\r),BC\r)==false disp('false');end
if isequal(double(bc\c),BC\c)==false disp('false');end
disp('end');
disp(' ');

%%{
BR\AR
br\ar

BR\AC
br\ac

BC\AR
bc\ar

BC\AC
bc\ac
%}

%{
1\AR
1\AC
1\ar
1\ac
%}

%{
BR\r
br\r

BR\c
br\c

BC\r
bc\r

BC\c
bc\c
%}

disp('=====================');
%{
A=[1 2 3];
A(:,:,2)=[4 5 6];
B=[2 4 6];
B(:,:,2)=[8 10 12];
A\B
%}

A=[1 2 3];
B=[2 4 6];
a=multi(A);
b=multi(B);
A\B
%a\b %この演算については、実装方法がわからないので、エラーメッセージを出力。

disp('=====================');

%%{
C=[1 2 3; 4 5 6; 7 8 9];
D=[10 11 12; 13 14 15; 16 17 18];
c=multi(C);
d=multi(D);
C\D
%c\d %この演算については、実装方法がわからないので、エラーメッセージを出力。

C=magic(3);
D=magic(3);
c=multi(C);
d=multi(D);
C\D
%c\d %この演算については、実装方法がわからないので、エラーメッセージを出力。
%}




