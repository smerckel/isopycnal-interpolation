#
#  $Id: Makefile,v 1e5e75c749c2 2015/08/08 22:03:51 fdelahoyde $
#
MAJOR_VERSION:=3
MINOR_VERSION:=05
PREFIX:=/usr/local

               CFLAGS:=	-O3
            CINCLUDES:=	-I.
              Library:=	libgswteos-10.so
              Program:=	gsw_check
      $(Program)_SRCS:=	gsw_check_functions.c \
			gsw_oceanographic_toolbox.c \
			gsw_saar.c

      $(Library)_SRCS:=	gsw_oceanographic_toolbox.c \
			gsw_saar.c 
      $(Library)_OBJS:=	gsw_oceanographic_toolbox.o \
			gsw_saar.o
      $(Library)_HDRS:=	gsw_internal_const.h gswteos-10.h

all:	$(Program)

$(Program):	$($(Program)_SRCS)
	gcc $(CFLAGS) $(CINCLUDES) -o $(Program) $($(Program)_SRCS) -lm

library:	$(Library)

$(Library):	$($(Library)_SRCS)
	gcc -fPIC -c $(CFLAGS) $(CINCLUDES) $($(Library)_SRCS)
	gcc -shared -Wl,-soname,$(Library).$(MAJOR_VERSION) -o $(Library).$(MAJOR_VERSION).$(MINOR_VERSION) $($(Library)_OBJS) -lm

clean:
	rm -f $(Program) $(Library) $($(Library)_OBJS)

install:
	install -m 0755 $(Library).$(MAJOR_VERSION).$(MINOR_VERSION) $(PREFIX)/lib
	ln -sf  $(PREFIX)/lib/$(Library).$(MAJOR_VERSION).$(MINOR_VERSION) $(PREFIX)/lib/$(Library)
	install -m 0755 $($(Library)_HDRS) $(PREFIX)/include
	ldconfig -n $(PREFIX)/lib


