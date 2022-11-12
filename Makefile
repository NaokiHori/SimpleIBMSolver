CC        := mpicc
CFLAGS    := -std=c99 -Wall -Wextra -O3 -flto -ffast-math -DNDIMS=2
DEPEND    := -MMD
LIBS      := -lfftw3 -lm
INCLUDES  := -Iinclude
SRCSDIR   := src
OBJSDIR   := obj
SRCS      := $(foreach dir, $(shell find $(SRCSDIR) -type d), $(wildcard $(dir)/*.c))
OBJS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.o)))
DEPS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.d)))
OUTPUTDIR := output
TARGET    := a.out

help:
	@echo "all      : create \"$(TARGET)\""
	@echo "clean    : remove \"$(TARGET)\" and object files \"$(OBJSDIR)/*.o\""
	@echo "output   : create directories to store artifacts"
	@echo "datadel  : remove contents in sub-directories of \"$(OUTPUTDIR)\""
	@echo "test_tdm : test tri-diagonam matrix solver"
	@echo "help     : show this help message"

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(DEPEND) -o $@ $^ $(LIBS)

$(OBJSDIR)/%.o: $(SRCSDIR)/%.c
	@if [ ! -e `dirname $@` ]; then \
		mkdir -p `dirname $@`; \
	fi
	$(CC) $(CFLAGS) $(DEPEND) $(INCLUDES) -c $< -o $@

clean:
	$(RM) -r $(OBJSDIR) $(TARGET)

output:
	@if [ ! -e $(OUTPUTDIR)/log ]; then \
		mkdir -p $(OUTPUTDIR)/log; \
	fi
	@if [ ! -e $(OUTPUTDIR)/save ]; then \
		mkdir -p $(OUTPUTDIR)/save; \
	fi
	@if [ ! -e $(OUTPUTDIR)/stat ]; then \
		mkdir -p $(OUTPUTDIR)/stat; \
	fi

datadel:
	$(RM) -r $(OUTPUTDIR)/log/*
	$(RM) -r $(OUTPUTDIR)/save/*
	$(RM) -r $(OUTPUTDIR)/stat/*

test_tdm:
	$(CC) -DDEBUG_TEST_TDM $(CFLAGS) $(INCLUDES) src/common.c src/tdm.c -o test_tdm.out -lm

-include $(DEPS)

.PHONY : help all clean output datadel test_tdm

