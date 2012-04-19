CXX = CC
CXXFLAGS = -O2

DIFF = ./sdiff
PRE = ./


.SUFFIXES:
.SUFFIXES: .a .o .c .cpp

.cpp.o:
		rm -f $*.cxx
		ln $*.cpp $*.cxx  
		$(CXX) $(CXXFLAGS) -c $*.cxx
		rm $*.cxx  

everything:    	tryrand 

tryrand_obj = tryrand.o newran.o myexcept.o extreal.o tryrand1.o tryrand2.o tryrand3.o tryrand4.o tryrand5.o hist.o

tryrand:       	$(tryrand_obj)
		$(CXX) -o $@ $(tryrand_obj) -L. -lm

tryrand.o:     	tryrand.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

newran.o:      	newran.cpp include.h newran.h boolean.h myexcept.h extreal.h

myexcept.o:    	myexcept.cpp include.h boolean.h myexcept.h

extreal.o:     	extreal.cpp include.h boolean.h extreal.h

tryrand1.o:    	tryrand1.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand2.o:    	tryrand2.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand3.o:    	tryrand3.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand4.o:    	tryrand4.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand5.o:    	tryrand5.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

hist.o:        	hist.cpp include.h boolean.h extreal.h newran.h tryrand.h myexcept.h

tryrand.txx:   	tryrand
		$(PRE)tryrand > tryrand.txx
		$(DIFF) tryrand.txt tryrand.txx

