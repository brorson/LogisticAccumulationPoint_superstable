CC := gcc
CPP := g++
CFLAGS := -w -Wtype-limits -Wextra
INCLUDEDIR := -I./eigen/ -I./mpfrc++
INCLUDES := globals.h

LIBLOCS := -L/usr/lib/x86_64-linux-gnu/
LDFLAGS := -lstdc++ -lmpfr -lgmp -lm -ldl


#SRCS := 
OBJS := superstable_calc.o secant.o secant_accelerated.o
EXES := superstable_calc

#=================================================

all: $(EXES)

# Compile
%.o: %.cpp $(INCLUDES)
	$(CPP) $(CFLAGS) $(INCLUDEDIR) -c $< -o $@

# Link 
superstable_calc: $(OBJS)
	echo "--> Linking ...."
	$(CPP) $(CFLAGS) $^ $(LIBLOCS) $(LDFLAGS) -o $@ 

# Clean
clean:
	-rm -f *.o *.obj *.out *.map $(EXES) $(OBJS) *~
