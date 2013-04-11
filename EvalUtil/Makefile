all: permtest

CXXFLAGS=-Wall -O5
#CXXFLAGS=-g


OBJS=permtest.o permtestmain.o compstat.o rand.o

$(OBJS) : *.hpp Makefile

permtest: $(OBJS)
	$(CXX) -o$@ $(OBJS) 

clean:
	rm -f *.o permtest

