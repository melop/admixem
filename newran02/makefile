#
# Makefile for 'newran'.
#
# Type 'make' or 'make newran' to create the binary.
# Type 'make clean' or 'make clear' to delete all temporaries.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# build target specs
CC = g++
CFLAGS = -O3 
OUT_DIR = release_build
LIBS =

# first target entry is the target invoked when typing 'make'
default: newran

newran: $(OUT_DIR)/extreal.cpp.o $(OUT_DIR)/hist.cpp.o $(OUT_DIR)/myexcept.cpp.o $(OUT_DIR)/newran.cpp.o
	@echo -n 'Linking newran... '
	@$(CC) $(CFLAGS) -o newran $(OUT_DIR)/extreal.cpp.o $(OUT_DIR)/hist.cpp.o $(OUT_DIR)/myexcept.cpp.o $(OUT_DIR)/newran.cpp.o $(LIBS)
	@echo Done.

$(OUT_DIR)/extreal.cpp.o: extreal.cpp include.h boolean.h extreal.h
	@echo -n 'Compiling extreal.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/extreal.cpp.o -c extreal.cpp
	@echo Done.

$(OUT_DIR)/hist.cpp.o: hist.cpp include.h boolean.h extreal.h newran.h \
 myexcept.h tryrand.h
	@echo -n 'Compiling hist.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/hist.cpp.o -c hist.cpp
	@echo Done.

$(OUT_DIR)/myexcept.cpp.o: myexcept.cpp include.h boolean.h myexcept.h
	@echo -n 'Compiling myexcept.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/myexcept.cpp.o -c myexcept.cpp
	@echo Done.

$(OUT_DIR)/newran.cpp.o: newran.cpp include.h newran.h boolean.h \
 myexcept.h extreal.h
	@echo -n 'Compiling newran.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/newran.cpp.o -c newran.cpp
	@echo Done.

run:
	./newran 

debug:
	gdb ./newran

clean:
	@echo -n 'Removing all temporary binaries... '
	@rm -f newran $(OUT_DIR)/*.o
	@echo Done.

clear:
	@echo -n 'Removing all temporary binaries... '
	@rm -f newran $(OUT_DIR)/*.o
	@echo Done.

