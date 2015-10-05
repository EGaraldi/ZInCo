#print some variables from each step of the program
#OPT += -DDEBUG

#print (a lot of) debug steps
#OPT += -DVDEBUG

#if your ICs use long IDs in the high-res ICs
#OPT += -DLONGIDS_IN_HIGH

#if your ICs use long IDs in the low-res ICs / snap
#OPT += -DLONGIDS_IN_LOW

#if you want to produce ICs with long IDs
#OPT += -DLONGIDS_OUT

#print (a lot of) more information to the info file (info file will become big!)
#OPT += -DVINFO

#produce a file whit the correspondence between old and new IDs
#OPT += -DPRINT_IDS_CORRESPONDENCE

CXX = mpic++

OPTIMIZE = -Wall -Wextra -g -O3

EXEC = ZInCo

OBJS = main.o  create_zoom_dilution_cascade.o  initialize_check_finalize_variables.o \
       New_files_writer.o  New_regions_finder.o  New_particles_maker.o \
       other_small_functions.o Slicer.o  Read_ICs_or_Snap_File.o  Read_parameter_file.o \
       Write_ICs_or_Snap_File.o  global_variables.o

INCL = function_declaration.h  global_variables.h  Makefile

CXXFLAGS = $(OPTIMIZE) $(OPT)

LIBS = -lm

%.o: %.cpp $(INCL)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EXEC): $(OBJS) 
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f $(OBJS) $(EXEC)
