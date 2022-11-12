CC        := cc
CFLAGS    := -O3 -std=c99 -flto -Wall -Wextra -DNDIMS=2
DEPEND    := -MMD
LIBS      := -lm
INCLUDES  := -Iinclude
SRCSDIR   := src
OBJSDIR   := obj
SRCS      := $(foreach dir, $(shell find $(SRCSDIR) -type d), $(wildcard $(dir)/*.c))
OBJS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.o)))
DEPS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.d)))
OUTPUTDIR := output
TARGET    := a.out

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files \"$(OBJSDIR)/*.o\""
	@echo "help    : show this help message"

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

-include $(DEPS)

.PHONY : help all clean

