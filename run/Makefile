include $(BASILISK)/Makefile.defs

CFLAGS+=-O2 -disable-dimensions $(BASILISK_INCLUDES) $(BASILISK_LIBS)

#HF1
slab-q1.c: slab-q.c
	ln -s slab-q.c slab-q1.c
slab-q1.s: CFLAGS += -Da_q=0.28 -Db_q=1.79 -Dtexp=665
slab-q1.tst:  CFLAGS += -Da_q=0.28 -Db_q=1.79 -Dtexp=665

#HF2
slab-q2.c: slab-q.c
	ln -s slab-q.c slab-q2.c
slab-q2.s: CFLAGS += -Da_q=0.44 -Db_q=1.79 -Dtexp=566
slab-q2.tst:  CFLAGS += -Da_q=0.44 -Db_q=1.79 -Dtexp=566

#HF3
slab-q3.c: slab-q.c
	ln -s slab-q.c slab-q3.c
slab-q3.s: CFLAGS += -Da_q=1.96 -Db_q=1.60 -Dtexp=479
slab-q3.tst:  CFLAGS += -Da_q=1.96 -Db_q=1.60 -Dtexp=479

#HF4
slab-q4.c: slab-q.c
	ln -s slab-q.c slab-q4.c
slab-q4.s: CFLAGS += -Da_q=1.13 -Db_q=1.72 -Dtexp=401
slab-q4.tst:  CFLAGS += -Da_q=1.13 -Db_q=1.72 -Dtexp=401

#HF5
slab-q5.c: slab-q.c
	ln -s slab-q.c slab-q5.c
slab-q5.s: CFLAGS += -Da_q=1.82 -Db_q=1.68 -Dtexp=360
slab-q5.tst:  CFLAGS += -Da_q=1.82 -Db_q=1.68 -Dtexp=360

