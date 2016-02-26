clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AR=[1 2 3;4 5 6;7 8 9];
AC=[1 2 3;4 5 6;7 8 9+i];

R=multi(AR);
C=multi(AC);

disp('start');

if isequal(double([R R]),[AR AR])==false disp('false');end
if isequal(double([R C]),[AR AC])==false disp('false');end
if isequal(double([C R]),[AC AR])==false disp('false');end
if isequal(double([C C]),[AC AC])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('start');

if isequal(double([C C C C R C C R]),[AC AC AC AC AR AC AC AR])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR=[1 2 3;4 5 6;7 8 9];
AC=[1 2;4 5;7 8+i];

R=multi(AR);
C=multi(AC);

disp('start');

if isequal(double([R R]),[AR AR])==false disp('false');end
if isequal(double([R C]),[AR AC])==false disp('false');end
if isequal(double([C R]),[AC AR])==false disp('false');end
if isequal(double([C C]),[AC AC])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('start');

if isequal(double([C C C C R C C R]),[AC AC AC AC AR AC AC AR])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR=[1 2 3;4 5 6;7 8 9];
AR(:,:,2)=[1 2 3;4 5 6;7 8 9];
AC=[1 2;4 5;7 8+i];
AC(:,:,2)=[1 2;4 5;7 8+i];

R=multi(AR);
C=multi(AC);

disp('start');

if isequal(double([R R]),[AR AR])==false disp('false');end
if isequal(double([R C]),[AR AC])==false disp('false');end
if isequal(double([C R]),[AC AR])==false disp('false');end
if isequal(double([C C]),[AC AC])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('start');

if isequal(double([C C C C R C C R]),[AC AC AC AC AR AC AC AR])==false disp('false');end

disp('end');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
