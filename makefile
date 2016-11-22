#written by Marvin Böttcher

SHELL= /bin/bash
CC:=g++ 



#compilation target name
TARGET := stochtreat

#directories for source code and build-files
BUILDDIR :=  build
SRCDIR :=  src
SRCEXT := cpp
INCLUDEDIR := include


#COMPILATION AND LIBARY FLAGS
# LDLIBS = -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas #-lm -lboost_serialization
LDLIBS = -L /usr/local/lib/ -lgsl #-lm -lboost_serialization

# -lefence -lboost_random -lboost_thread-mt  -lboost_random  -stdlib=libc++

WARNINGS := -Wall -pedantic -Wextra -Wshadow  -Wcast-qual  -Weffc++ -Wfloat-equal -Wunreachable-code   -Wdisabled-optimization 
#--Wstrict-overflow=5 Wpointer-arith -Wunused -Wconversion -fno-diagnostics-fixit-info

OPTIDEBUG := -O3
#CFLAGS: release mode: -O2 or -O3 debug mode: -g -ggdb  profile mode: -pg 

CFLAGS := -std=c++11 $(OPTIDEBUG) $(WARNINGS) -I./$(INCLUDEDIR)/


#formatting of variables##################################################
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEPS = $(OBJECTS:%.o=%.d)

# $(info make sources: $(SOURCES))
# $(info make objects: $(OBJECTS))
# echo $(SOURCES)

#OBJSUT = $(patsubst ./build/main.o, ,$(OBJS)) $(addprefix $(BUILDDIR),$(OBJUnitTest))
#DEPSUT = $(OBJSUNITTEST:%.o=%.d)

#building all objects (main program) #########################################  
all: $(TARGET)
# rm tags
# ctags -R .
	
#build the main target
$(TARGET): $(OBJECTS)
	$(CC)  $(CFLAGS)  -o $(TARGET) $(OBJECTS) $(LDLIBS)


#all .o rules  ################################################### 
./$(BUILDDIR)/%.o : ./$(SRCDIR)/%.$(SRCEXT)
	$(CC) $(CFLAGS)  -c -o $@ $<



#copy files to "/data/"
copy: all
	cp $(TARGET) ./data/


#create dependencies ############################################
dep: $(DEPS)

-include $(DEPS) #include depencies (saved in .d files)

#rule for creating the dependency files###############################################
$(BUILDDIR)/%.d: $(SRCDIR)/%.$(SRCEXT)
	$(CC) $(CFLAGS) -MM  -MT '$(subst .d,.o,$@)'  $< -MF $@ 



#clean ##################################################################
.PHONY : clean

clean : 
	rm -f $(TARGET) $(OBJECTS) $(DEPS) $(OBJSUT) $(OBJSPT)


#run the compiled program#############################
.PHONY: run

run: all
	./$(TARGET)


