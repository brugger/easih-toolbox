# En Makefile som den burde være lavet fra starten
# Lavet til luna og venner.
# Lavet af Kim Brugger Sommeren 2000
PROG   = tile2seq
WARN   = #-Wshadow  -Wmissing-prototypes
FLAGS  = -pipe -Wall -ggdb -I$(ILIB) -D_POSIX_SOURCE $(WARN) -O6  #-pg
CC     = gcc
LIBS   = #-lm #-lefence #-ldmalloc  #-lefence #-lmcheck#-lmalloc
RM     = rm -f
MAKE   = gmake
SRCS   = tile2seq
OBJS   = tile2seq.o


all :  $(PROG)  

tile2seq:  tile2seq.o 
	$(CC) $(FLAGS)  -o  $@  $<   $(LIBS)


%.o : %.c
	$(CC) $(FLAGS) -c  $< -o $@

clean: 
	$(RM) $(PROG)

veryclean : cleandoc
	$(RM) $(ILIB)*.o *.o core $(PROG) *~ $(ILIB)*~ *.bak\
# DO NOT DELETE
