#!/usr/bin/env python
#/*=========================================================================
#
#    Copyright (c) 2008-2011 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#    See License_v1.2.txt for details.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/
#
# Script to patch the ParaView source tree for builtin reader installation.
#

import os, shutil, subprocess, sys

if len(sys.argv) != 2:
    print 'Usage: ', sys.argv[0], ' <top directory of ParaView source tree>'
    sys.exit(1)

ParaViewPath = sys.argv[1] + '/'
vtkPath = ParaViewPath + 'VTK/'

print 'Copying files'

vtkIOPath = vtkPath + 'IO'
shutil.copy2("../vtkOFFDevReader.cxx", vtkIOPath)
shutil.copy2("../vtkOFFDevReader.h", vtkIOPath)

vtkParallelPath = vtkPath + 'Parallel'
shutil.copy2("../vtkPOFFDevReader.cxx", vtkParallelPath)
shutil.copy2("../vtkPOFFDevReader.h", vtkParallelPath)

vtkCMakePath = vtkPath + 'CMake'
shutil.copy2("../CMake/FindRegex.cmake", vtkCMakePath)

pqComponentsPath = ParaViewPath + 'Qt/Components'
shutil.copy2("../pqPOFFDevReaderPanel.cxx", pqComponentsPath)
shutil.copy2("../pqPOFFDevReaderPanel.h", pqComponentsPath)

print 'Editing XML'

resourcesPath = ParaViewPath + 'Servers/ServerManager/Resources/'
readersXMLFile = resourcesPath + 'readers.xml'
readersXMLHead, readersXMLTail = open(readersXMLFile).read() \
    .split('    <!-- Beginning of OpenFOAM Reader -->')
dummy, readersXMLTail = readersXMLTail \
    .split('    <!-- End of OpenFOAM Reader -->')

smXMLHead, smXMLContents = open('../OFFReaderSM.xml').read() \
    .split('  <ProxyGroup name="sources">')
smXMLContents, dummy = smXMLContents.split('  </ProxyGroup>')

readersXMLHandle = open(readersXMLFile, 'wb')
readersXMLHandle.write(readersXMLHead)
readersXMLHandle.write(smXMLContents)
readersXMLHandle.write(readersXMLTail)
readersXMLHandle.close()

print 'Applying patch'

patchHandle = open('CMakeLists.diff')
os.chdir(ParaViewPath)
patchProcess = subprocess.Popen(['patch', '-p1'], stdin=patchHandle)
patchProcess.wait()
patchHandle.close()

print 'Done.'
