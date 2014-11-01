C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = MGL.o

 MGL.so: $(OBJECTS)
	R CMD SHLIB -o MGL.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so
