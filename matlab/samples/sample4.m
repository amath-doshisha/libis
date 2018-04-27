%function sample4
clear all

%% initializations
tic
N=2;
M=4;
disp('【パラメータの設置】');
disp(sprintf('帯幅: N=%d M=%d',N,M));
lambda=[1 2 3 4 5 6 7 8 9 10];
disp(sprintf('指定固有値: [%s]',num2str(lambda)));
%prec=53;lambda=double(lambda);
prec=512;set_default_prec(prec);lambda=multi(lambda);
disp(sprintf('計算精度: %d bits',prec));

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
A=sample_matTN(lambda,N,M);
disp('生成された行列: A=');
disp(num2str(A,'%.2f     '));
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp('倍精度にキャストされた行列double(A)の固有値を倍精度のeigで計算: ');
l=sort(eig(double(A)));
disp(num2str(l,'%.40f'))
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp(sprintf('生成された行列Aを多倍長(%d bits)のeigで計算:',prec*64));
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp(sprintf('倍精度にキャストされた行列double(A)を多倍長(%d bits)のeigで計算:',prec*64));
A=double(A);
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
toc