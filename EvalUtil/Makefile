all: permtest

# gcc/clang optimization flags
#CXXFLAGS+=-Wall -Ofast
# gcc optimization flags
CXXFLAGS+=-Wall -O3
# add debug info if required
#CXXFLAGS+=-g


OBJS=permtest.o permtestmain.o compstat.o rand.o

$(OBJS) : *.hpp Makefile

permtest: $(OBJS)
	$(CXX) -o$@ $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o permtest

