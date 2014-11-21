# CMake script for SubStory package
# Author: Bastien Boussau
# Created: 21/11/2014

ADD_EXECUTABLE(substory SubStory.cpp)
TARGET_LINK_LIBRARIES(substory ${LIBS})

# Install progs
INSTALL(TARGETS
  substory
  DESTINATION bin)

