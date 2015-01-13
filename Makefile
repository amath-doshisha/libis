
all:
	( cd src ; make )
	( cd test ; make )
	( cd prog ; make )
	( cd matlab ; make )
	( cd docs ; make )
clean:
	( cd src ; make clean )
	( cd test ; make clean )
	( cd prog ; make clean )
	( cd matlab ; make clean )
	( cd docs ; make clean )
	rm -fv *~
