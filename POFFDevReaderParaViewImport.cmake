#/*=========================================================================
#
#    Copyright (c) 2009-2010 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

INCLUDE(${POFFDevReader_SOURCE_DIR}/CMake/FindRegex.cmake)

SET(POFFDevReader_SRCS 
  ${POFFDevReader_SOURCE_DIR}/vtkOFFDevReader.cxx
  ${POFFDevReader_SOURCE_DIR}/vtkPOFFDevReader.cxx
  )

INCLUDE_DIRECTORIES(${POFFDevReader_SOURCE_DIR})
#INCLUDE_DIRECTORIES(${POFFDevReader_SOURCE_DIR}/..)

PARAVIEW_INCLUDE_WRAPPED_SOURCES("${POFFDevReader_SRCS}")

PARAVIEW_INCLUDE_SERVERMANAGER_RESOURCES(
  "${POFFDevReader_SOURCE_DIR}/POFFDevReaderSM.xml"
  )
