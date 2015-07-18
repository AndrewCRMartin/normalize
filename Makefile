HOME = /home/bsm/martin
COPT = -I$(HOME)/include
LOPT = -L$(HOME)/lib
CC = cc -ansi -pedantic -Wall
OFILES1 = normalize.o erf.o
OFILES2 = gendata.o
TIFILES = algorithm.aux algorithm.dvi algorithm.log
all : normalize gendata algorithm.pdf


normalize : $(OFILES1)
	$(CC) $(LOPT) -o $@ $(OFILES1) -lgen -lm

gendata : $(OFILES2)
	$(CC) $(LOPT) -o $@ $(OFILES2) -lgen -lm

algorithm.pdf : algorithm.tex
	latex algorithm
	dvipdf algorithm

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f $(OFILES1) $(OFILES2) $(TIFILES)