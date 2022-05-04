COMPILER := gfortran
CFLAGS := -O1

PRG := main
SRC := $(wildcard *.f90 *.f)
OBJ := $(patsubst %.f, %.o, $(patsubst %.f90, %.o, $(SRC)))

build: $(PRG)

$(PRG): $(OBJ)
	$(COMPILER) $(OBJ) -o $@

main.o: modd.o dop853.o

%.o: %.f90
	$(COMPILER) $(CFLAGS) -c -o $@ $<

%.o: %.f
	$(COMPILER) $(CFLAGS) -c -o $@ $<

res: $(PRG)
	./$<
	spd-say 'toots toots toots'
	
plot: main
	python3 plotting.py

clean:
	rm -f *.o *.mod $(PRG) f_RESULT E_RESULT
	echo "Cleaned!"

.PHONY: clean res build plot
.SILENT: res clean
