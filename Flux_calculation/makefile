
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
         -I $(ROOT_INC) \
         -I $(HEP_CONCURRENCY_INC) \
         -I $(LARSIM_DIR)


CXXFLAGS=-std=c++17 -Wall
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
	-L build_slf7.x86_64/larsim/lib -l larsim_EventWeight_Base


EXEC=plot


$(EXEC): joseph/GeoVector.cxx joseph/GeoAABox.cxx joseph/GeoHalfLine.cxx joseph/GeoLine.cxx joseph/GeoLineSegment.cxx joseph/GeoCone.cxx joseph/GeoSphere.cxx joseph/GeoTrajectory.cxx joseph/GeoAlgo.cxx plot.c 
	 @echo Building $(EXEC)
	 @$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
   
#$(EXEC): plot.c
#	 @echo Building $(EXEC)
#	 @$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
