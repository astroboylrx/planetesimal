CC = g++

CFLAGS = -c -Wall -g

OPTIMIZE = -O2

LDFLAGS = -g

LIBS   = #-lm -lgsl -lgslcblas

SOURCES = disruption.cpp main.cpp coagulation.cpp fragmentation.cpp condensation.cpp cog_frag_ingrav.cpp pre_define.cpp

OBJECTS = $(SOURCES:.cpp=.o)

EXEC = planetesimal

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CC) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS)  $(LIBS) -o $@

.cpp.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $< -o $@

.PHONY : clean

clean:
	 rm -f $(OBJECTS) $(EXEC)
