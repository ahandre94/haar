CC = g++
CFLAGS = -O2 -Wall -fopenmp
LDFLAGS = -lgomp -lpthread `pkg-config --cflags --libs opencv` -lstdc++fs

TARGET = haar

SRCDIR = src
OBJDIR = obj
IMGDIR = output_img
IMGDIRSEQ = output_img_seq
FSTATS = stats.txt

SOURCES = $(shell find $(SRCDIR) -type f -name *.cpp)
SRCDIRS = $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJ = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOURCES))
HEADER = $(shell find $(SRCDIR) -type f -name *.h)

all: directories $(TARGET)

directories:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(IMGDIR)
	@mkdir -p $(IMGDIRSEQ)

$(OBJDIR)/$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(TARGET): buildrepo $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

$(OBJ): $(HEADER)

buildrepo:
	@$(call make-repo)

define make-repo
	for dir in $(SRCDIRS); \
	do \
		mkdir -p $(OBJDIR)/$$dir; \
	done
endef

clean:
	rm -rf $(OBJDIR)

cleanall:
	rm -rf $(TARGET) $(OBJDIR) $(IMGDIR) $(IMGDIRSEQ) $(FSTATS)

.PHONY:
	all clean cleanall
