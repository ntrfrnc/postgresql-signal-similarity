PG_CONFIG = pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
INCLUDEDIR = $(shell $(PG_CONFIG) --includedir-server)
INSTALLDIR = $(shell $(PG_CONFIG) --pkglibdir)
include $(PGXS)

all: signal_similarity.so

signal_similarity.so: signal_similarity.o
	cc -O3 -shared -o signal_similarity.so signal_similarity.o

signal_similarity.o: signal_similarity.c
	cc -O3 -o signal_similarity.o -c signal_similarity.c $(CFLAGS) -I$(INCLUDEDIR)

install: signal_similarity.so
	install -m 644 signal_similarity.so $(INSTALLDIR)/

uninstall: 
	rm -f $(INSTALLDIR)/signal_similarity.so

clean:
	rm -f signal_similarity.o signal_similarity.so
