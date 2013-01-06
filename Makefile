include Makefile.amd64-ubuntu

DEST		= .

LINKER		= $(ARCHLINKER)

LDFLAGS 	= $(ARCHLDFLAGS) 

LIBS    	= $(ARCHLIBS)

DEBUG_CFLAGS	= $(ARCHDEBUGFLAGS) $(COMMON_CFLAGS) -DDEBUG

PROFILE_CFLAGS	= $(ARCHPROFFLAGS)

OPENGL_CFLAGS	= $(DEBUG_CFLAGS) -DUSE_OPENGL_VIEWER

COMMON_CFLAGS	= -Inumtools

CFLAGS		= $(ARCHCFLAGS) $(COMMON_CFLAGS)

CCFLAGS		= $(ARCHCCFLAGS)

CC		= $(ARCHCC)
CCC		= $(ARCHCCC)

SRCS		= main.c propagate.c queue.c read_comp_react.c\
		  stats.c analyse.c random.c


OBJS		= $(ARCHOBJS) \
	          main.o propagate.o queue.o read_comp_react.o\
		  stats.o analyse.o random.o

OUTPUT		= Gillespie
OUTPUT2		= Gillespie.2

vers2:  $(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -o $(OUTPUT2)

all:	$(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -o $(OUTPUT)

opengl:
	$(MAKE) "CFLAGS=$(OPENGL_CFLAGS)" all

debug:
	$(MAKE) "CFLAGS=$(DEBUG_CFLAGS)" all

profile:
	$(MAKE) "CFLAGS=-g -pg" pall

pall:	$(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -g -pg -o $(OUTPUT)

depend: 
	makedepend -f Makefile.powerpc -- $(CFLAGS) -- $(SRCS)

clean:;	@echo "Cleaning... "
	@rm -f *.o "#"* core *.bak
	@echo "Done."



# DO NOT DELETE
