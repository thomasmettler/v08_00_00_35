# Master build script for CMake (automatically generated by "mrb newDev")

# Require minimum cmake version
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

# Add top level include
include_directories ( ${PROJECT_SOURCE_DIR} )

# Enable testing
enable_testing()

# add_subdirectory and other commands will be added below by @mrb gitCheckout@ for each package

# If you add a package by hand, you MUST also add the appropriate package block below
# Remove the appropriate package block below if you delete a product from your srcs area
# Or, use mrb uc to update this file safely

# DO NOT DELETE 
# ubcrt package block
set(ubcrt_not_in_ups true)
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/ubcrt )
include_directories ( $ENV{MRB_BUILDDIR}/ubcrt )

# DO NOT DELETE 
ADD_SUBDIRECTORY(ubcrt)

