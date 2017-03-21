CC = g++
CFLAGS =-O3 -g -std=c++11
DEPS = simulation.hpp interaction.hpp polymer.hpp point.hpp parameters.hpp observable.hpp numgen.hpp unsymmeig.hpp symmeig.hpp svd.hpp GLE-OPT.hpp timer.hpp bias.hpp spline.hpp
OBJS = main.o simulation.o interaction.o polymer.o point.o parameters.o observable.o numgen.o unsymmeig.o symmeig.o svd.o GLE-OPT.o bias.o spline.o

pimd: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

%.o : %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o
