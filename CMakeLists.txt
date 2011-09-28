#/*=========================================================================
#
#    Copyright (c) 2008-2010 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#    See License_v1.2.txt for details.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Find a ParaView build first in order to get its configuration

IF(ParaView_SOURCE_DIR)
  SET(INSOURCE_PLUGIN TRUE)
ELSE()
  SET(INSOURCE_PLUGIN FALSE)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF()

# Construct version number string and check

IF(NOT DEFINED PARAVIEW_VERSION)
  SET(PARAVIEW_VERSION "${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}")
ENDIF()

IF(PARAVIEW_VERSION VERSION_LESS "3.8")
  MESSAGE(FATAL_ERROR "ParaView version ${PARAVIEW_VERSION} is too old."
    " The package only supports ParaView 3.8 or later.")
ENDIF()

IF(NOT INSOURCE_PLUGIN)

  # Set default installation directory

  IF(PARAVIEW_LIBRARY_DIRS)
    IF(CMAKE_INSTALL_PREFIX STREQUAL "/usr/local"
        OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files/Project"
        OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files (x86)/Project")
      IF(DEFINED ENV{PV_PLUGIN_PATH})
        # Installing to plugin path defined by PV_PLUGIN_PATH
        SET(CMAKE_INSTALL_PREFIX "$ENV{PV_PLUGIN_PATH}")
      ELSE()
        IF(ParaView_BINARY_DIR MATCHES ".*/bin")
          # Installing to the development install of ParaView
          IF(APPLE)
            SET(CMAKE_INSTALL_PREFIX "${ParaView_BINARY_DIR}/Plugins")
          ELSEIF(DEFINED CMAKE_BUILD_TYPE)
            SET(CMAKE_INSTALL_PREFIX "${ParaView_BINARY_DIR}/plugins")
          ELSE()
            # This is not checked yet; may be wrong
            SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_BINARY_DIR}/Release/plugins")
          ENDIF()
        ELSE()
          # Installing to the build directory of ParaView
          IF(APPLE)
            SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/Plugins")
          ELSEIF(DEFINED CMAKE_BUILD_TYPE)
            SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/plugins")
          ELSE()
            SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/Release/plugins")
          ENDIF()
        ENDIF()
      ENDIF()
      SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
        "Install path prefix, prepended onto install directories." FORCE)
    ENDIF()
  ENDIF()


  # Set default build type to Release if Makefile-based build system is used

  IF(DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  ENDIF()

ENDIF()

# THE SERVER SIDE #############################################################

# Add include directories

IF(INSOURCE_PLUGIN)
  INCLUDE_DIRECTORIES(${PARAVIEW_INCLUDE_DIRS} ${VTK_INCLUDE_DIR})
ENDIF()

# Find a regular expression library

#SET(Project_CMAKE_DIR "${Project_SOURCE_DIR}/CMake")
#INCLUDE(${Project_CMAKE_DIR}/FindRegex.cmake)
INCLUDE("CMake/FindRegex.cmake")

# For ParaView 3.9 or greater, disable locale workaround since it's
# been fixed in ParaView itself

IF(PARAVIEW_VERSION VERSION_GREATER "3.8")
  ADD_DEFINITIONS(-DVTK_FOAMFILE_LOCALE_WORKAROUND=0)
ENDIF()

# Server side plugin

ADD_PARAVIEW_PLUGIN(
  POFFDevReaderPlugin "1.0"
  REQUIRED_ON_SERVER
  SERVER_MANAGER_XML OFFReaderSM.xml
  SERVER_MANAGER_SOURCES vtkOFFDevReader.cxx vtkPOFFDevReader.cxx
  )

# Set installation destination

IF(NOT INSOURCE_PLUGIN)
  # ParaView 3.6.2 has been modified after release so that plugins are
  # built without manifests if MSVC is used. With the modification
  # there seems to be no way to tell exactly whether plugins are built
  # as SHARED or MODULE so there seems to be no way to tell whether
  # RUNTIME or LIBRARY target should be specified either.
  INSTALL(TARGETS POFFDevReaderPlugin DESTINATION "${CMAKE_INSTALL_PREFIX}")
ENDIF()

# THE CLIENT SIDE #############################################################

IF(PARAVIEW_BUILD_QT_GUI)

  # Add include directories

  IF(INSOURCE_PLUGIN)
    INCLUDE_DIRECTORIES(
      ${PARAVIEW_GUI_INCLUDE_DIRS}
    )
    INCLUDE(${QT_USE_FILE})
  ENDIF()

  # Additional UI elements for the OFFReader panel

  QT4_WRAP_CPP(
    MOC_SRCS pqPOFFDevReaderPanel.h
    )

  ADD_PARAVIEW_OBJECT_PANEL(
    IFACES IFACE_SRCS
    CLASS_NAME pqPOFFDevReaderPanel
    XML_NAME OpenFOAMReader
    XML_GROUP sources
    )

  ADD_PARAVIEW_PLUGIN(
    POFFDevReaderPanelPlugin "1.0"
    REQUIRED_ON_CLIENT
    GUI_INTERFACES ${IFACES}
    GUI_SOURCES pqPOFFDevReaderPanel.cxx
    ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS}
    )

  # Set installation destination

  IF(NOT INSOURCE_PLUGIN)
    # ParaView 3.6.2 has been modified after release so that plugins
    # are built without manifests if MSVC is used. With the
    # modification there seems to be no way to tell exactly whether
    # plugins are built as SHARED or MODULE so there seems to be no
    # way to tell whether RUNTIME or LIBRARY target should be
    # specified either.
    INSTALL(TARGETS POFFDevReaderPanelPlugin DESTINATION
      "${CMAKE_INSTALL_PREFIX}")
  ENDIF()
ENDIF()
