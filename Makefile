INC_DIR=-I/home/anl/include/ -I/home/stafik/gl/anl/include/
LIB_DIR=-L/home/anl/lib/ -L/home/stafik/gl/anl/lib/

CFLAGS=-std=c++11 -no-pie -g

default: main.out

%.o: %.cpp
	g++ $(CFLAGS) -c $(INC_DIR)  $<

%.out: %.o
	g++ $(CFLAGS) $<  -o  $@   -lGLEW -lGL -lglfw   -L$(LIB_DIR) -lcommon

clean:
	rm -f *.o *.out *~
