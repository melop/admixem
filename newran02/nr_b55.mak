BORLANDPATH = "C:\program files\Borland\cbuilder5"

TASM = TASM32
TLIB = tlib
TLINK = ilink32
LIBPATH = $(BORLANDPATH)\LIB
INCLUDEPATH = $(BORLANDPATH)\INCLUDE

DIFF = sdiff
PRE =

CC = bcc32 -W- -v- -H- -3 -N -Og -Oi -Ov -f -I$(INCLUDEPATH)

.cpp.obj:
   $(CC) -c {$< }

everything:    	tryrand.exe 

tryrand_obj = tryrand.obj newran.obj myexcept.obj extreal.obj tryrand1.obj tryrand2.obj tryrand3.obj tryrand4.obj tryrand5.obj hist.obj

tryrand.exe:   	$(tryrand_obj)
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(tryrand_obj),$@,,import32.lib cw32.lib
|

tryrand.obj:   	tryrand.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

newran.obj:    	newran.cpp include.h newran.h boolean.h myexcept.h extreal.h

myexcept.obj:  	myexcept.cpp include.h boolean.h myexcept.h

extreal.obj:   	extreal.cpp include.h boolean.h extreal.h

tryrand1.obj:  	tryrand1.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand2.obj:  	tryrand2.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand3.obj:  	tryrand3.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand4.obj:  	tryrand4.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

tryrand5.obj:  	tryrand5.cpp include.h newran.h tryrand.h boolean.h myexcept.h extreal.h

hist.obj:      	hist.cpp include.h boolean.h extreal.h newran.h tryrand.h myexcept.h

tryrand.txx:   	tryrand.exe
		$(PRE)tryrand > tryrand.txx
		$(DIFF) tryrand.txt tryrand.txx

