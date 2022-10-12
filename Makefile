.SUFFICES: .C

#CC = g++ -g -Wall -pedantic -Wno-undefined-var-template
CC = g++ -O3 -Wall -pedantic -funroll-loops -Wno-undefined-var-template
#CC = g++ -O3 -Wall -pedantic -funroll-loops 
#CC = icpc  -O3 -funroll-loops -static

# Please set the SPRNG_ROOT environment variable to the
# directory for the SPRNG libraries.
# INC_SPRNG = -I$(SPRNG_ROOT)/SRC
# LIB_SPRNG = -L$(SPRNG_ROOT) -L$(SPRNG_ROOT)/lib -lsprng

# INCLUDES = $(INC_SPRNG)
#LIBS = $(LIB_SPRNG) -lm -static
LIBS = $(LIB_SPRNG) -lm
FLAGS = -c

#Compiler rule for c++ compiler.
.C.o :
	$(CC) $(FLAGS) $(INCLUDES) $<
.cc.o :
	$(CC) $(FLAGS) $(INCLUDES) $<

CMC_OBJ =  Random.o Standard.o PeriodicPoint.o BasicPoint.o AccRatio.o\
	CMC.o Sampling.o PrintErrorBar.o InitialPos.o Histogram.o Physics.o\
	Timer.o Spline.o Form.o ReadInTable.o Parser.o LSpline.o\
	State.o PrintDiffusion.o WriteVASPPOSCARFile.o

cmc: $(CMC_OBJ)
	$(CC) $(CMC_OBJ) $(LIBS) -o $@


clean:	
	rm -f *.o *.ii *.ti cmc cci dal h2int

depend :
	perl ./listdepend.scr dependencies.txt

include dependencies.txt
