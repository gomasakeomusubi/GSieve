CC=g++
OPTFLAGS= -march=native -O2 -pipe
STRICT= -std=c++98 -pedantic -Wall -fno-builtin
CFLAGS   = $(STRICT) $(OPTFLAGS)
LIBS = -lm -lgmp -lntl -lpthread

OBJ_DIR = obj
BIN_DIR = bin
BASE = tool sampler gsieve
MAIN  = main

##
TEMP_BASE_OBJ = $(addsuffix .o, $(BASE)  )
BASE_OBJ = $(addprefix $(OBJ_DIR)/, $(TEMP_BASE_OBJ) )

all: $(MAIN)

$(OBJ_DIR)/%.o: %.cc
	$(CC) $(CONFOPT) $(CFLAGS) -c   $< -o $@

$(MAIN): $(BASE_OBJ)
	$(CC) $@.cc $(BASE_OBJ) -o $(BIN_DIR)/$@  $(LIBS) $(CFLAGS)

clean:
	rm bin/* obj/*