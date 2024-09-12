TSK_VER=1.1.2
TARBALL=tskit-${TSK_VER}.tar.xz
DOWNLOAD=https://github.com/tskit-dev/tskit/releases/download/C_${TSK_VER}/${TARBALL}

TSK_OBJECTS=kastore.o tskit_tables.o tskit_core.o tskit_trees.o \
	tskit_stats.o tskit_genotypes.o tskit_convert.o
CFLAGS = -Wall -Wextra -O3 -Itskit/subprojects/kastore -Itskit
LDFLAGS= -lgsl

run: run.c tree.c recorder.c mcmc.c libtskit.a
	${CC} ${CFLAGS} $^ -o $@ -lm ${LDFLAGS}
	
verification: verification.c tree.c recorder.c mcmc.c libtskit.a
	${CC} ${CFLAGS} $^ -o $@ -lm ${LDFLAGS}

libtskit.a: ${TSK_OBJECTS}
	${AR} rcs $@ ${TSK_OBJECTS}

kastore.o: tskit
	${CC} -c ${CFLAGS} tskit/subprojects/kastore/kastore.c -o kastore.o

tskit_%.o: tskit
	${CC} -c ${CFLAGS} tskit/tskit/$*.c -o $@

${TARBALL}:
	wget ${DOWNLOAD}

tskit: ${TARBALL}
	tar xf ${TARBALL}
	mv tskit-${TSK_VER} tskit

clean:
	rm -f *.a *.o

mrproper: clean
	rm -fr tskit example example.trees tskit-*.tar.xz
