CC=gfortran
CFLAGS = -O3 -ffree-form -ffast-math 
all: 
	$(CC) $(CFLAGS) -c xmulti_common.f
	$(CC) $(CFLAGS) -c parse_text.f
	$(CC) $(CFLAGS) -c allocate.f
	$(CC) $(CFLAGS) -c nr_common.f
	$(CC) $(CFLAGS) -c nr_mod.f
	$(CC) $(CFLAGS) xmulti.f xmulti_common.o nr_mod.o nr_common.o allocate.o parse_text.o -o xmulti
	mv ./xmulti ../bin/.
