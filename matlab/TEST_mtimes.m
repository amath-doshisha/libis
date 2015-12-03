clear all
clc


%A=[1 2 3;4 5 6;7 8 9];
%A(:,:,2)=[10 20 30;40 50 60;70 80 90];

AR=[1 2 3;4 5 6;7 8 9];
AC=[1 2 3;4 5 6;7 8 9+i];

R=multi(AR);
C=multi(AC);

ar=1;
ac=1+i;

r=multi(ar);
c=multi(ac);

disp('start');

if isequal(double(R*R),AR*AR)==false disp('false');end
if isequal(double(R*C),AR*AC)==false disp('false');end
if isequal(double(C*R),AC*AR)==false disp('false');end
if isequal(double(C*C),AC*AC)==false disp('false');end

if isequal(double(r*R),ar*AR)==false disp('false');end
if isequal(double(r*C),ar*AC)==false disp('false');end
if isequal(double(c*R),ac*AR)==false disp('false');end
if isequal(double(c*C),ac*AC)==false disp('false');end

if isequal(double(R*r),AR*ar)==false disp('false');end
if isequal(double(R*c),AR*ac)==false disp('false');end
if isequal(double(C*r),AC*ar)==false disp('false');end
if isequal(double(C*c),AC*ac)==false disp('false');end

disp('end');
disp(' ');




clear all
disp('start');

AR1=[1 2 3;4 5 6];
AR2=[1 4;2 5;3 6];

AC1=[1 2 3;4 5 6+i];
AC2=[1 4;2 5;3 6+i];

R1=multi(AR1);
R2=multi(AR2);

C1=multi(AC1);
C2=multi(AC2);

if isequal(double(R1*R2),AR1*AR2)==false disp('false');end
if isequal(double(R2*R1),AR2*AR1)==false disp('false');end

if isequal(double(R1*C2),AR1*AC2)==false disp('false');end
if isequal(double(R1*C2),AR1*AC2)==false disp('false');end
disp('end');


clear all
















