clear all
set_default_prec(53)
x=multi(665857);
y=multi(470832);
disp('The variables x and y are storaged by 53-bit.');
disp('x='); disp(x);
disp('y='); disp(y);
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 53-bit.');
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is NOT true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 64-bit.');
set_default_prec(64)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is NOT true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 128-bit.');
set_default_prec(128)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 256-bit.');
set_default_prec(256)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 53-bit and automatic precision.');
set_default_prec(53);
auto_prec_enabled();
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is NOT true.');
disp('[Important] Automatic precision works for only operators "+", "-" and "*".');
disp('------');
disp('Compute z=x*x*x*x-4*y*y-4*y*y*y*y by 53-bit and automatic precision.');
z=x*x*x*x-4*y*y-4*y*y*y*y;
disp('z='); disp(z);
disp('This result is true.');
disp(sprintf('The precison of z is %d-bit.',get_prec(z)));
disp('------');
disp('Compute the following values by 53-bit and automatic precision.');
z1=x*x*x*x;
z2=-4*y*y;
z3=-4*y*y*y*y;
z=z1+z2+z3;
disp(sprintf('z1=x*x*x*x=%s',num2str(z1,'%.0f')));
disp(sprintf('z2=-4*y*y=%s',num2str(z2,'%.0f')));
disp(sprintf('z3=-4*y*y*y*y=%s',num2str(z3,'%.0f')));
disp(sprintf('z=z1+z2+z3=%s',num2str(z,'%.0f')));
disp('These results are true.');
disp(sprintf('The precison of z1 is %d-bit.',get_prec(z1)));
disp(sprintf('The precison of z2 is %d-bit.',get_prec(z2)));
disp(sprintf('The precison of z3 is %d-bit.',get_prec(z3)));
disp(sprintf('The precison of z is %d-bit.',get_prec(z)));

auto_prec_disabled()