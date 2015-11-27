clc

AR=[1 2 3;4 5 6;7 8 9];
AC=[1 2 3;4 5 6;7 8 9+i];

ar=1;
ac=1+i;

r=multi(ar);
c=multi(ac);

R=multi(AR);
C=multi(AC);

disp('start');

if isequal(double(R-R),AR-AR)==false disp('false');end
if isequal(double(R-C),AR-AC)==false disp('false');end
if isequal(double(C-R),AC-AR)==false disp('false');end
if isequal(double(C-C),AC-AC)==false disp('false');end

if isequal(double(r-R),ar-AR)==false disp('false');end
if isequal(double(r-C),ar-AC)==false disp('false');end
if isequal(double(c-R),ac-AR)==false disp('false');end
if isequal(double(c-C),ac-AC)==false disp('false');end

if isequal(double(R-r),AR-ar)==false disp('false');end
if isequal(double(R-c),AR-ac)==false disp('false');end
if isequal(double(C-r),AC-ar)==false disp('false');end
if isequal(double(C-c),AC-ac)==false disp('false');end

disp('end');

%{
R-R
R-C
C-R
C-C

r-R
r-C
c-R
c-C

R-r
R-c
C-r
C-c
%}















