# this is my make file
Target = Genes

SOURCES = main.cpp SNP_Fun.cpp

OBJECTS= $(SOURCES:.cpp=.o)

CXXFLAGS = $-std=c++11 -Wall -Wextra -I

all: ${TARGET}

clean:
		rm ${OBJECTS} ${TARGET}

${TARGET}: ${OBJECTS}
		${LINK.cc} -o $@ $^
