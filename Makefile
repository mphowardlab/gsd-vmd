# set this to: /path/to/your/vmd
VMDDIR=/Applications/VMD\ 1.9.3.app/Contents/vmd

CFLAGS=-m32 -O3 -Wall -pedantic -std=c99 -fPIC -I$(VMDDIR)/plugins/include
LDFLAGS=-dynamiclib -L/usr/X11R6/lib

.PHONY: all clean
all: gsdplugin.so

TARGET=gsdplugin.so
SOURCES=gsdplugin.c gsd.c
OBJECTS=$(SOURCES:.c=.o)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS) -lm

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
