function y=imag(x)
y=rmulti(-1,cmulti_data(get_prec,'one','imag',x.data));
