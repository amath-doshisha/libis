clear all
set_default_prec(53)
x=multi(665857);
y=multi(470832);

disp('[53-bit]');
z=x.^4-4*y.^2-4*y.^4
disp('[64-bit]');
set_default_prec(64)
z=x.^4-4*y.^2-4*y.^4
disp('[128-bit]');
set_default_prec(128)
z=x.^4-4*y.^2-4*y.^4
disp('[256-bit]');
set_default_prec(256)
z=x.^4-4*y.^2-4*y.^4

disp('[auto precision mode]');
set_default_prec(53)
auto_prec_enabled()
z=x*x*x*x-4*y*y-4*y*y*y*y
z.get_prec

y=multi(470831);
z=x*x*x*x-4*y*y-4*y*y*y*y
z.get_prec
