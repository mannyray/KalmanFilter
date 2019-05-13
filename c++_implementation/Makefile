# text followed by '#' is commented out
CC := g++ 
SRCDIR := src
BUILDDIR := build
TARGET := bin/runner
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -std=c++11 -O2  # -Wall
LIB := -L lib  
INC := -I include -I /usr/local/include/eigen3/

all: setup sample_code 

setup:
	$(shell mkdir -p bin)

sample_code:
	$(shell mkdir -p bin/sample_code)
	$(shell cp ../testing/testing_continuousDiscreteKalmanFilter/*.txt bin/sample_code)
	$(CC) $(CFLAGS) test/sample_code/sample_code.cc $(SOURCES) $(INC) $(LIB) -o bin/sample_code/sample_code

sample_code2:
	$(shell mkdir -p bin/sample_code2)
	$(shell cp ../testing/testing_discreteDiscreteKalmanFilter/*.txt bin/sample_code2)
	$(CC) $(CFLAGS) test/sample_code2/sample_code2.cc $(SOURCES) $(INC) $(LIB) -o bin/sample_code2/sample_code2


$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
