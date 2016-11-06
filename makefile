all:
	make clean; make program; make run; make plot

# dla C++ wszedzie zmieniamy .c na .cpp

CC= g++		# dla C++:   CC=g++
CFLAGS=
INCLUDE=
LIB= -lm	# dla lapacka:	LIB= -lm -llapack -lblas
SOURCES= 
OBJECTS= $(SOURCES:.cpp=.o)

TARGET = program

$(TARGET): main.cpp $(OBJECTS)
	$(CC) -o $@ $^ $(LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

run:
	./program

plot:
	gnuplot 'plot.gp'


.PHONY: clean

clean:
	rm -f *.o *.txt *.png $(TARGET)

