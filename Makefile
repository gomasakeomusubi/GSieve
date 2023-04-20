CC            = g++
CFLAGS        = 
DEST          = 
LDFLAGS       = -lntl -lgmp -lpthread -lgsl -lgslcblas 
LIBS          = 
OBJS          = main.o
PROGRAM       = main
SRC = main.cpp
HEADER = sample.h gsieve.h

all:$(PROGRAM)

$(PROGRAM):$(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

$(OBJS):$(SRC) $(HEADER)
	$(CC) $(SRC) $(LDFLAGS) $(LIBS) -c

# get_basis:get_basis.cpp
# 	$(CC) get_basis.cpp -lcurl -o get_basis

clean:;         
	rm -f *.o *~ $(PROGRAM)

install:$(PROGRAM)
	install -s $(PROGRAM) $(DEST)
retags:;
	rm TAGS
	etags -a *.cpp *.h