SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES:.cpp .o
OBJS = btt.o
PROG = btt
CC = g++
DEBUG = -O3
#DEBUG = -g -D DEBUG -O0
OPT = -flto  
CFLAGS = -Wall $(OPT)  $(DEBUG)
LDFLAGS = -Wall  $(OPT) $(DEBUG)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(PROG)

btt.o: btt.cpp
	$(CC) -I. -I$(srcdir)$(CFLAGS) -c $< -o $@

clean:
	$(RM) $(PROG) $(OBJS)
