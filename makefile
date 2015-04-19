#
# Makefile for 'admixem'.
#
# Type 'make' or 'make admixem' to create the binary.
# Type 'make clean' or 'make clear' to delete all temporaries.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# build target specs
CC = g++
CFLAGS = -fopenmp -O3 
OUT_DIR = release_build
NEW_RAN_DIR = ./newran02
PARSER_DIR = ./parser
LIBS =
MKDIR_P = mkdir -p


# first target entry is the target invoked when typing 'make'
default: directories admixemp

directories: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}



admixemp: $(OUT_DIR)/parser.cpp.o $(OUT_DIR)/extreal.cpp.o $(OUT_DIR)/hist.cpp.o $(OUT_DIR)/myexcept.cpp.o $(OUT_DIR)/newran.cpp.o $(OUT_DIR)/config.cpp.o $(OUT_DIR)/filesystem.cpp.o $(OUT_DIR)/Individual.cpp.o $(OUT_DIR)/main.cpp.o $(OUT_DIR)/makemarkerfile.cpp.o $(OUT_DIR)/maths.cpp.o $(OUT_DIR)/Population.cpp.o $(OUT_DIR)/simulation.cpp.o $(OUT_DIR)/testexpression.cpp.o
	@echo -n 'Linking admixem... '
	@$(CC) $(CFLAGS) -o bin/admixemp $(OUT_DIR)/parser.cpp.o $(OUT_DIR)/extreal.cpp.o $(OUT_DIR)/hist.cpp.o $(OUT_DIR)/myexcept.cpp.o $(OUT_DIR)/newran.cpp.o  $(OUT_DIR)/config.cpp.o $(OUT_DIR)/filesystem.cpp.o $(OUT_DIR)/Individual.cpp.o $(OUT_DIR)/main.cpp.o $(OUT_DIR)/makemarkerfile.cpp.o $(OUT_DIR)/maths.cpp.o $(OUT_DIR)/Population.cpp.o $(OUT_DIR)/simulation.cpp.o $(OUT_DIR)/testexpression.cpp.o $(LIBS)
	@cp scripts/*.php bin/
	@echo Done.
	
	

$(OUT_DIR)/config.cpp.o: config.cpp config.h parser/parser.h \
 newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h maths.h makemarkerfile.h
	@echo -n 'Compiling config.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/config.cpp.o -c config.cpp
	@echo Done.

$(OUT_DIR)/filesystem.cpp.o: filesystem.cpp
	@echo -n 'Compiling filesystem.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/filesystem.cpp.o -c filesystem.cpp
	@echo Done.

$(OUT_DIR)/Individual.cpp.o: Individual.cpp Individual.h makemarkerfile.h \
 newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h config.h parser/parser.h maths.h
	@echo -n 'Compiling Individual.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/Individual.cpp.o -c Individual.cpp
	@echo Done.

$(OUT_DIR)/main.cpp.o: main.cpp makemarkerfile.h testexpression.h \
 simulation.h config.h parser/parser.h newran02/newran.h \
 newran02/include.h newran02/boolean.h newran02/myexcept.h \
 newran02/extreal.h maths.h Population.h Individual.h
	@echo -n 'Compiling main.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/main.cpp.o -c main.cpp
	@echo Done.

$(OUT_DIR)/makemarkerfile.cpp.o: makemarkerfile.cpp def.h \
 makemarkerfile.h newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h maths.h
	@echo -n 'Compiling makemarkerfile.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/makemarkerfile.cpp.o -c makemarkerfile.cpp
	@echo Done.

$(OUT_DIR)/maths.cpp.o: maths.cpp maths.h makemarkerfile.h \
 newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h
	@echo -n 'Compiling maths.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/maths.cpp.o -c maths.cpp
	@echo Done.

$(OUT_DIR)/Population.cpp.o: Population.cpp Population.h Individual.h \
 makemarkerfile.h newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h config.h parser/parser.h maths.h
	@echo -n 'Compiling Population.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/Population.cpp.o -c Population.cpp
	@echo Done.

$(OUT_DIR)/simulation.cpp.o: simulation.cpp simulation.h config.h \
 parser/parser.h newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h maths.h makemarkerfile.h \
 Population.h Individual.h
	@echo -n 'Compiling simulation.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/simulation.cpp.o -c simulation.cpp
	@echo Done.

$(OUT_DIR)/testexpression.cpp.o: testexpression.cpp testexpression.h \
 parser/parser.h
	@echo -n 'Compiling testexpression.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/testexpression.cpp.o -c testexpression.cpp
	@echo Done.

	
$(OUT_DIR)/extreal.cpp.o: $(NEW_RAN_DIR)/extreal.cpp $(NEW_RAN_DIR)/include.h $(NEW_RAN_DIR)/boolean.h $(NEW_RAN_DIR)/extreal.h
	@echo -n 'Compiling extreal.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/extreal.cpp.o -c $(NEW_RAN_DIR)/extreal.cpp
	@echo Done.

$(OUT_DIR)/hist.cpp.o: $(NEW_RAN_DIR)/hist.cpp $(NEW_RAN_DIR)/include.h $(NEW_RAN_DIR)/boolean.h $(NEW_RAN_DIR)/extreal.h $(NEW_RAN_DIR)/newran.h \
 $(NEW_RAN_DIR)/myexcept.h $(NEW_RAN_DIR)/tryrand.h
	@echo -n 'Compiling hist.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/hist.cpp.o -c $(NEW_RAN_DIR)/hist.cpp
	@echo Done.

$(OUT_DIR)/myexcept.cpp.o: $(NEW_RAN_DIR)/myexcept.cpp $(NEW_RAN_DIR)/include.h $(NEW_RAN_DIR)/boolean.h $(NEW_RAN_DIR)/myexcept.h
	@echo -n 'Compiling myexcept.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/myexcept.cpp.o -c $(NEW_RAN_DIR)/myexcept.cpp
	@echo Done.

$(OUT_DIR)/newran.cpp.o: $(NEW_RAN_DIR)/newran.cpp $(NEW_RAN_DIR)/include.h $(NEW_RAN_DIR)/newran.h $(NEW_RAN_DIR)/boolean.h \
 $(NEW_RAN_DIR)/myexcept.h $(NEW_RAN_DIR)/extreal.h
	@echo -n 'Compiling newran.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/newran.cpp.o -c $(NEW_RAN_DIR)/newran.cpp
	@echo Done.
	
$(OUT_DIR)/parser.cpp.o: $(PARSER_DIR)/parser.cpp $(PARSER_DIR)/parser.h
	@echo -n 'Compiling parser.cpp... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/parser.cpp.o -c $(PARSER_DIR)/parser.cpp
	@echo Done.

run:
	./admixemp 

debug:
	gdb ./admixemp

clean:
	@echo -n 'Removing all temporary binaries... '
	@rm -f admixemp $(OUT_DIR)/*.o
	@echo Done.

clear:
	@echo -n 'Removing all temporary binaries... '
	@rm -f admixemp $(OUT_DIR)/*.o
	@echo Done.

