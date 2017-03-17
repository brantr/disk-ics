EXEC   = disk_ics

OPTIMIZE =  -O2  

OBJS   = main.o disk_ics.o rk_int.o

CC     = g++

INCL   = disk_ics.h rk_int.h

LIBS   = -lm

CFLAGS = $(OPTIMIZE)

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)   

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

