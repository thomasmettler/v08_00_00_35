#Completely stolen.
CPPFLAGS=-I $(BOOST_INC) \
	 -I $(CANVAS_ROOT_IO_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
	 -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(NUSIMDATA_INC) \
	 -I $(LARSIM_INC) \
         -I $(ROOT_INC)

CXXFLAGS=-std=c++14 -Wall
CXX=g++
LDFLAGS=$$(root-config --libs) \
	-L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io\
	-L $(CANVAS_LIB) -l canvas\
        -L $(CETLIB_LIB) -l cetlib \
	-L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
     -L ../build_slf7.x86_64/larsim/lib -l larsim_EventWeight_Base
	###-L $(LARSIM_LIB) -l larsim_EventWeight_Base


UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)
  EXEC=example_program
else
  EXEC=example_program
endif

$(EXEC): GeoVector.cxx GeoAABox.cxx GeoHalfLine.cxx GeoLine.cxx GeoLineSegment.cxx GeoCone.cxx GeoSphere.cxx GeoTrajectory.cxx GeoAlgo.cxx test.cc 
	 @echo Building $(EXEC)
	 @$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
