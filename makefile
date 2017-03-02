CC = g++
CFLAGS =-O3 -g -std=c++11
DEPS = simulation.hpp interaction.hpp polymer.hpp point.hpp parameters.hpp observable.hpp
OBJS = main.o simulation.o interaction.o polymer.o point.o parameters.o observable.o

pimd: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

%.o : %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o
