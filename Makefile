CC := gcc
CPP := g++
CFLAGS := -w
INCLUDEDIR := -I./eigen/ -I./mpfrc++
INCLUDES := f.h J.h globals.h

LIBLOCS := -L/usr/lib/x86_64-linux-gnu/
LDFLAGS := -lstdc++ -lmpfr -lgmp


#SRCS := 
OBJS := bifurcation_calc.o  f.o J.o
EXES := bifurcation_calc

#=================================================

all: $(EXES)

%.o: %.cpp $(INCLUDES)
	$(CC) $(CFLAGS) $(INCLUDEDIR) -c $< -o $@

# Link 
bifurcation_calc: $(OBJS)
	echo "--> Linking ...."
	$(CC) $(CFLAGS) $^ $(LIBLOCS) $(LDFLAGS) -o $@ 


# Clean
clean:
	-rm -f *.o *.obj *.out *.map $(EXES) $(OBJS) *~
