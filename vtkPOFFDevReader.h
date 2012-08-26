/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPOFFReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See License_v1.2.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPOFFReader - reads a decomposed dataset in OpenFOAM(R) file format
// .SECTION Description
// vtkPOFFReader creates a multiblock dataset. It reads
// parallel-decomposed mesh information and time dependent data.  The
// polyMesh folders contain mesh information. The time folders contain
// transient data for the cells. Each folder can contain any number of
// data files.

// .SECTION Thanks
// This class was developed by Takuya Oshima at Niigata University,
// Japan (oshima@eng.niigata-u.ac.jp).

// .SECTION Disclaimer
// OPENFOAM(R) is a registered trade mark of OpenCFD Limited, the
// producer of the OpenFOAM software and owner of the OPENFOAM(R) and
// OpenCFD(R) trade marks. This code is not approved or endorsed by
// OpenCFD Limited.



#ifndef __vtkPOFFDevReader_h
#define __vtkPOFFDevReader_h

#include "vtkOFFDevReader.h"

class vtkDataArray;
class vtkDataArraySelection;
class vtkDoubleArray;
class vtkIntArray;
class vtkMultiProcessController;
class vtkPolyData;
class vtkStringArray;
class vtkUnstructuredGrid;

class
#if !defined(POFFDevReaderPlugin_EXPORTS)
VTK_PARALLEL_EXPORT
#endif
vtkPOFFDevReader : public vtkOFFDevReader
{
public:
  //BTX
  enum caseType
    {
    DECOMPOSED_CASE_APPENDED = 0,
    RECONSTRUCTED_CASE = 1,
    DECOMPOSED_CASE_MULTIBLOCK = 2
    };
  //ETX
  static vtkPOFFDevReader *New();

  vtkTypeMacro(vtkPOFFDevReader, vtkOFFReader);

  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Set and get case type. 0 = decomposed case (appended), 1 = reconstructed
  // case, 2 = decomposed case (multiblock).
  void SetCaseType(const int);
  vtkGetMacro(CaseType, caseType);

  void SetCaseTypeToDecomposedCaseAppended()
  { this->SetCaseType(DECOMPOSED_CASE_APPENDED); }

  void SetCaseTypeToReconstructedCase()
  { this->SetCaseType(RECONSTRUCTED_CASE); }

  void SetCaseTypeToDecomposedCaseMultiblock()
  { this->SetCaseType(DECOMPOSED_CASE_MULTIBLOCK); }

  // Description:
  // Get enum value that represents reconstructed case.
  static caseType GetReconstructedCase() { return RECONSTRUCTED_CASE; }

  // Description:
  // Get enum value that represents decomposed case.
  static caseType GetDecomposedCaseAppended()
  { return DECOMPOSED_CASE_APPENDED; }

  // Description:
  // Get enum value that represents reconstructed case (multiblock).
  static caseType GetDecomposedCaseMultiblock()
  { return DECOMPOSED_CASE_MULTIBLOCK; }

  // Description:
  // Set and Get if region (patch) names should be displayed in the
  // client render window. Currenly only works for reconstructed cases.
  vtkSetMacro(ShowRegionNames, int);
  vtkGetMacro(ShowRegionNames, int);
  vtkBooleanMacro(ShowRegionNames, int);

  // Description:
  // Set and get the communication object used to exchange metadata
  // between the processes.
  virtual void SetController(vtkMultiProcessController *);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  // Description:
  // Get if region names should be redrawn.
  vtkGetMacro(RedrawRegionNames, int);

  // Description:
  // Get the array that holds region centroids.
  vtkGetObjectMacro(RegionCentroids, vtkDoubleArray);

  // Description:
  // Get the array that holds region names.
  vtkGetObjectMacro(RegionNames, vtkStringArray);

  // Description:
  // Get the array that holds region name display styles.
  vtkGetObjectMacro(RegionStyles, vtkIntArray);

  // Description:
  // Instruct the reader that the updated region names has been drawn by the
  // client side.
  void RedrewRegionNames()
  {
    this->RedrawRegionNames = 0;
  }

protected:
  vtkPOFFReader();
  ~vtkPOFFReader();

  int RequestInformation(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);

private:
  vtkMultiProcessController *Controller;
  caseType CaseType;
  int MaximumNumberOfPieces;
  int NumProcesses;
  int ProcessId;
  vtkIntArray *PieceIds;
  int MaximumPieceId; // must be signed

  int ShowRegionNames;
  int RedrawRegionNames;
  vtkDoubleArray *RegionCentroids;
  vtkStringArray *RegionNames;
  vtkIntArray *RegionStyles;

  caseType CaseTypeOld;
  unsigned long MTimeOld;
  int ShowRegionNamesOld;
  int NReadersOld;

  vtkPOFFReader(const vtkPOFFReader &); // Not implemented.
  void operator=(const vtkPOFFReader &); // Not implemented.

  void CreateRegionArrays(vtkMultiBlockDataSet *);
  void MultiBlockRegionCentroids(vtkMultiBlockDataSet *, const char *);
  int InternalRegionCentroid(vtkUnstructuredGrid *);
  int AppendedPolyRegionCentroids(vtkPolyData *);
  int PolyRegionCentroids(vtkPolyData *);

  void GatherMetaData();
  void AllReduceStatus(int &);
  void Broadcast(vtkStringArray *, int);
  void ConstructBlocks(vtkMultiBlockDataSet *, const int *, int &,
      vtkStringArray *, vtkIdType &);
  void BroadcastStructure(vtkMultiBlockDataSet *, const int);
  void GatherV(vtkStringArray *);
  void AllGatherV(vtkStringArray *, const bool);
  void AllGatherV(vtkDataArraySelection *);
  void GatherV(vtkDataArray *, vtkDataArray *);
  void AllGatherV(vtkDataArray *, vtkDataArray *);
};

#endif
