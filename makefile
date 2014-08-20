EXEC   = planetesimal

OPTIMIZE  =  -O2  -m64 

OPTIMIZE += $(OPTS)

OBJS  =  main.o coagulation.o fragmentation.o condensation.o cog_frag_ingrav.o pre_define.o

CC     = mpicxx

INCL   = coagulation.h fragmentation.h condensation.h cog_frag_ingrav.h pre_define.h

LIBS   = -lmpicxx -lmpi -lpmpi -lpthread

CFLAGS  = $(OPTIMIZE) -I/opt/local/include/mpich-mp/
CPPFLAGS  = $(OPTIMIZE) -I/opt/local/include/mpich-mp/
LDFLAGS = -m64 

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LDFLAGS)  $(LIBS) -o $(EXEC)
	 rm -f $(OBJS)   

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)
     

