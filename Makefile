CC = gcc
CFLAGS = -Wall -pedantic -Werror
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	main.o\
      	visual.o\
      	logger.o\
      	boundary_configurator.o


all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(OBJ)

helper.o      : helper.h logger.h
init.o        : helper.h init.h boundary_configurator.h logger.h
boundary_val.o: helper.h boundary_val.h logger.h
uvp.o         : helper.h uvp.h logger.h
visual.o      : helper.h logger.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h logger.h boundary_configurator.h
