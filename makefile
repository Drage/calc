COMPILER = g++
PROGRAM = calc
OPTIONS = -O2

# Search for code files
SRC += $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SRC))

all : $(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(COMPILER) $(OPTIONS) -o $(PROGRAM) $(OBJS)
	rm *.o

%.o : %.c
	$(COMPILER) -c $<
