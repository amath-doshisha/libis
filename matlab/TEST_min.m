clear all
clc

a=[1 2 3;4 5 6;7 8 9];
a(:,:,2)=[10 20 30;40 50 60;70 80 90+i];
A=multi(a)

disp('=================');

min(A)
min(min(A))
min(min(min(A)))


angle([-1-i,1-i,1+i,-1+i])%•¡‘f”‚Åâ‘Î’l‚ª“™‚µ‚¢ê‡‚ÍAˆÊ‘ŠŠp‚Ì¬‚³‚¢‡(‚±‚ê‚Ì¬‚³‚¢‡)
min([-1-i,1-i,1+i,-1+i])
min(multi([-1-i,1-i,1+i,-1+i]))
min(multi([-1-i,1-i,1+i]))
min(multi([-1-i,1-i]))
