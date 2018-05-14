
all:
	( cd src ; ./compile.sh )
	( cd test ; make clean ; make )
#	( cd prog ; make )
	( cd matlab ; make )

clean:
	rm -rfv ./MACOSX ./LINUX ./include
	rm -fv *~
	( cd src ; make clean )
	( cd test ; make clean )
	( cd prog ; make clean )
	( cd matlab ; make clean )
	( cd docs ; make clean )
