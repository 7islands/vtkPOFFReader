/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPOFFReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See License_v1.2.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This class was developed by Takuya Oshima at Niigata University,
// Japan (oshima@eng.niigata-u.ac.jp).
// Disclaimer:
// OPENFOAM(R) is a registered trade mark of OpenCFD Limited, the
// producer of the OpenFOAM software and owner of the OPENFOAM(R) and
// OpenCFD(R) trade marks. This code is not approved or endorsed by
// OpenCFD Limited.

#include "vtkPOFFDevReader.h"

#include "vtkAppendCompositeDataLeaves.h"
#include "vtkAppendPolyData.h"
#include "vtkCharArray.h"
#include "vtkCleanPolyData.h"
#include "vtkCollection.h"
#include "vtkDataArraySelection.h"
#include "vtkDataObjectTypes.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkSortDataArray.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVariantArray.h"
#include "vtkVersion.h"

#if VTK_MAJOR_VERSION >= 6
#include "vtkDataObjectTreeIterator.h"
#else
#include "vtkCompositeDataIterator.h"
#endif

#include <vtksys/ios/sstream>

vtkStandardNewMacro(vtkPOFFReader);
vtkCxxSetObjectMacro(vtkPOFFReader, Controller, vtkMultiProcessController);

//-----------------------------------------------------------------------------
vtkPOFFReader::vtkPOFFReader()
{
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL)
    {
    this->NumProcesses = 1;
    this->ProcessId = 0;
    }
  else
    {
    this->NumProcesses = this->Controller->GetNumberOfProcesses();
    this->ProcessId = this->Controller->GetLocalProcessId();
    }
  this->CaseType = RECONSTRUCTED_CASE;
  this->MTimeOld = 0;
  this->MaximumNumberOfPieces = 0;
  this->PieceIds = vtkIntArray::New();
  this->MaximumPieceId = -1;
  this->NReadersOld = -1;

  this->ShowRegionNames = 0;
  this->ShowRegionNamesOld = 0;
  this->RedrawRegionNames = 0;
  this->RegionCentroids = vtkDoubleArray::New();
  this->RegionCentroids->SetNumberOfComponents(3);
  this->RegionNames = vtkStringArray::New();
  this->RegionStyles = vtkIntArray::New();

  this->CaseTypeOld = RECONSTRUCTED_CASE;
}

//-----------------------------------------------------------------------------
vtkPOFFReader::~vtkPOFFReader()
{
  this->PieceIds->Delete();
  this->RegionCentroids->Delete();
  this->RegionNames->Delete();
  this->RegionStyles->Delete();

  this->SetController(NULL);
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Controller: " << this->Controller << endl;
  os << indent << "Case Type: " << this->CaseType << endl;
  os << indent << "Show Region Names: " << this->ShowRegionNames << endl;
  os << indent << "Redraw Region Names: " << this->RedrawRegionNames << endl;
  os << indent << "Maximum Number of Pieces: " << this->MaximumNumberOfPieces
      << endl;
  os << indent << "Number of Processes: " << this->NumProcesses << endl;
  os << indent << "Process Id: " << this->ProcessId << endl;
  os << indent << "Piece Ids: \n";
  this->PieceIds->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Maximum Piece Id: " << this->MaximumPieceId << endl;
  os << indent << "MTimeOld: " << this->MTimeOld << endl;
  os << indent << "CaseTypeOld: " << this->CaseTypeOld << endl;
  os << indent << "Show Region Names Old: " << this->ShowRegionNamesOld << endl;
  os << indent << "NReaders Old: " << this->NReadersOld << endl;
  os << indent << "RegionCentroids: \n";
  this->RegionCentroids->PrintSelf(os, indent.GetNextIndent());
  os << indent << "RegionNames: \n";
  this->RegionNames->PrintSelf(os, indent.GetNextIndent());
  os << indent << "RegionStyles: \n";
  this->RegionStyles->PrintSelf(os, indent.GetNextIndent());
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::SetCaseType(const int t)
{
  if (this->CaseType != t)
    {
    this->CaseType = static_cast<caseType>(t);
    this->Refresh = true;
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
int vtkPOFFReader::RequestInformation(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  this->SetErrorCode(vtkErrorCode::NoError);
  int ret = 1;

  if (this->CaseType == RECONSTRUCTED_CASE)
    {
    // refresh selection lists when switching from DECOMPOSED_CASE_MULTIBLOCK
    if (this->CaseTypeOld == DECOMPOSED_CASE_MULTIBLOCK)
      {
      // clear patch selections
      this->Superclass::PatchDataArraySelection->RemoveAllArrays();
      }

    this->CaseTypeOld = this->CaseType;
    this->Superclass::SetOutputProcessorPatches(PROCESSOR_PATCHES_ON);

    if (this->ProcessId == 0)
      {
      ret = this->Superclass::RequestInformation(request, inputVector,
          outputVector);
      }
    if (this->NumProcesses > 1)
      {
      // if there was an error in process 0 abort all processes
      this->AllReduceStatus(ret);
      if (ret == 0)
        {
        vtkErrorMacro(<< "The master process returned an error.");
        this->SetErrorCode(vtkErrorCode::UnknownError);
        return 0;
        }

      vtkDoubleArray *timeValues;
      if (this->ProcessId == 0)
        {
        timeValues = this->Superclass::GetTimeValues();
        }
      else
        {
        timeValues = vtkDoubleArray::New();
        }
      this->Controller->Broadcast(timeValues, 0);
      if (this->ProcessId != 0)
        {
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        timeValues->Delete();
        this->Superclass::Refresh = false;
        }
      this->Controller->Broadcast(this->CasePath, 0);
      this->GatherMetaData(); // pvserver deadlocks without this
      }

    return ret;
    }

  if (!this->Superclass::FileName || strlen(this->Superclass::FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    this->SetErrorCode(vtkErrorCode::NoFileNameError);
    return 0;
    }

  if (*this->Superclass::FileNameOld != this->Superclass::FileName
      || this->Superclass::ListTimeStepsByControlDict
          != this->Superclass::ListTimeStepsByControlDictOld
      || this->Superclass::Refresh)
    {
    // retain selection statuses when just refreshing a case
    if (*this->Superclass::FileNameOld != "" && *this->Superclass::FileNameOld != this->Superclass::FileName)
      {
      // clear selections
      this->Superclass::CellDataArraySelection->RemoveAllArrays();
      this->Superclass::SurfaceDataArraySelection->RemoveAllArrays();
      this->Superclass::PointDataArraySelection->RemoveAllArrays();
      this->Superclass::LagrangianDataArraySelection->RemoveAllArrays();
      this->Superclass::PatchDataArraySelection->RemoveAllArrays();
      }
    else if ((this->CaseTypeOld == DECOMPOSED_CASE_MULTIBLOCK
            && this->CaseType != DECOMPOSED_CASE_MULTIBLOCK)
        || (this->CaseTypeOld != DECOMPOSED_CASE_MULTIBLOCK
        && this->CaseType == DECOMPOSED_CASE_MULTIBLOCK))
      {
      this->Superclass::PatchDataArraySelection->RemoveAllArrays();
      }

    // this substitution is always run because Refresh is always set
    // when the CaseType changes
    this->CaseTypeOld = this->CaseType;

    *this->Superclass::FileNameOld = vtkStdString(this->FileName);
    this->Superclass::Readers->RemoveAllItems();
    this->Superclass::NumberOfReaders = 0;
    this->PieceIds->Initialize();
    this->MaximumPieceId = -1;
    this->MaximumNumberOfPieces = 0;
    this->Superclass::SetOutputProcessorPatches(
        this->CaseType == DECOMPOSED_CASE_MULTIBLOCK
        ? PROCESSOR_PATCHES_AGGREGATE : PROCESSOR_PATCHES_OFF);

    vtkIntArray *procNos = vtkIntArray::New();
    vtkStringArray *procNames = vtkStringArray::New();
    vtkDoubleArray *timeValues;

    // recreate case information
    vtkStdString masterCasePath, controlDictPath;
    this->Superclass::CreateCasePath(masterCasePath, controlDictPath);

    this->Superclass::CreateCharArrayFromString(this->Superclass::CasePath,
        "CasePath", masterCasePath);

    if (this->ProcessId == 0)
      {
      // search and list processor subdirectories
      vtkDirectory *dir = vtkDirectory::New();
      if (!dir->Open(masterCasePath.c_str()))
        {
        vtkErrorMacro(<< "Can't open " << masterCasePath.c_str());
        this->SetErrorCode(vtkErrorCode::UnknownError);
        procNos->Delete();
        procNames->Delete();
        dir->Delete();
        this->AllReduceStatus(ret = 0);
        return 0;
        }
      for (int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
        {
        const vtkStdString subDir(dir->GetFile(fileI));
        if (subDir.substr(0, 9) == "processor")
          {
          const vtkStdString procNoStr(subDir.substr(9, vtkStdString::npos));
          char *conversionEnd;
          const int procNo = strtol(procNoStr.c_str(), &conversionEnd, 10);
          if (procNoStr.c_str() + procNoStr.length() == conversionEnd && procNo
              >= 0)
            {
            procNos->InsertNextValue(procNo);
            procNames->InsertNextValue(subDir);
            }
          }
        }
      procNos->Squeeze();
      procNames->Squeeze();
      dir->Delete();

      // sort processor subdirectories by processor numbers
      vtkSortDataArray::Sort(procNos, procNames);

      // The reader does not allow processor directories with
      // duplicated processor number (e.g. processor1 and
      // processor01). On the contrary missing processor subdirectory
      // is allowed.
      for (int procI = 1; procI < procNos->GetNumberOfTuples(); procI++)
        {
        if (procNos->GetValue(procI - 1) == procNos->GetValue(procI))
          {
          vtkWarningMacro(<< "Processor subdirectories with duplicated "
              "processor number " << procNames->GetValue(procI - 1).c_str()
              << " and " << procNames->GetValue(procI).c_str() << " found. "
              << procNames->GetValue(procI).c_str() << " will be ignored.");
          procNos->RemoveTuple(procI);
          // vtkStringArray does not have RemoveTuple()
          for (int procJ = procI + 1; procJ < procNames->GetNumberOfTuples();
              procJ++)
            {
            procNames->SetValue(procJ - 1, procNames->GetValue(procJ));
            }
          procNames->Resize(procNames->GetNumberOfTuples() - 1);
          }
        }

      // get time directories from the first processor subdirectory
      if (procNames->GetNumberOfTuples() > 0)
        {
        vtkOFFReader *masterReader = vtkOFFReader::New();
        masterReader->SetFileName(this->FileName);
        masterReader->SetParent(this);
        if (!masterReader->MakeInformationVector(outputVector, procNames
            ->GetValue(0)) || !masterReader->MakeMetaDataAtTimeStep(true))
          {
          procNos->Delete();
          procNames->Delete();
          masterReader->Delete();
          this->AllReduceStatus(ret = 0);
          return 0;
          }
        this->Superclass::Readers->AddItem(masterReader);
        this->PieceIds->InsertNextValue(procNos->GetValue(0));
        timeValues = masterReader->GetTimeValues();
        masterReader->Delete();
        }
      else
        {
        timeValues = vtkDoubleArray::New();
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        }
      }
    else
      {
      timeValues = vtkDoubleArray::New();
      }

    if (this->NumProcesses > 1)
      {
      // if there was an error in process 0 abort all processes
      this->AllReduceStatus(ret);
      if (ret == 0)
        {
        vtkErrorMacro(<< "The master process returned an error.");
        this->SetErrorCode(vtkErrorCode::UnknownError);
        timeValues->Delete(); // don't have to care about process 0
        procNos->Delete();
        procNames->Delete();
        return 0;
        }

      this->Controller->Broadcast(procNos, 0);
      this->Broadcast(procNames, 1);
      this->Controller->Broadcast(timeValues, 0);
      if (this->ProcessId != 0)
        {
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        timeValues->Delete();
        }
      }

    if (this->ProcessId == 0 && procNames->GetNumberOfTuples() == 0)
      {
      vtkWarningMacro(<< "Case " << this->CasePath->GetPointer(0) << " contains no processor subdirectory.");
      timeValues->Delete();
      }

    this->MaximumNumberOfPieces = procNames->GetNumberOfTuples();

    if (procNos->GetNumberOfTuples() > 0)
      {
      this->MaximumPieceId
          = procNos->GetValue(procNos->GetNumberOfTuples() - 1);
      }
    else
      {
      this->MaximumPieceId = -1;
      }

    // create reader instances for other processor subdirectories
    // skip processor0 since it's already been created
    for (int procI = (this->ProcessId ? this->ProcessId : this->NumProcesses); procI
        < procNames->GetNumberOfTuples(); procI += this->NumProcesses)
      {
      vtkOFFReader *subReader = vtkOFFReader::New();
      subReader->SetFileName(this->FileName);
      subReader->SetParent(this);
      ret *= (subReader->MakeInformationVector(NULL, procNames->GetValue(procI))
          && subReader->MakeMetaDataAtTimeStep(true));
      // subreader is added even if it failed in order to match the
      // max number of pieces and the total number of readers
      this->Superclass::Readers->AddItem(subReader);
      // does not use procI but use processor subdirectory number so
      // that the block number corresponds to processor subdirectory
      // number. By this way missing processor subdirectories can
      // easily be identified.
      this->PieceIds->InsertNextValue(procNos->GetValue(procI));
      subReader->Delete();
      }

    procNos->Delete();
    procNames->Delete();

    this->PieceIds->Squeeze();

    this->GatherMetaData();
    this->Superclass::Refresh = false;
    }

  // it seems MAXIMUM_NUMBER_OF_PIECES must be returned every time
  // RequestInformation() is called
  outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
      this->MaximumNumberOfPieces);
  return ret;
}

//-----------------------------------------------------------------------------
int vtkPOFFReader::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet
      *output =
          vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int ret = 1;

  if (this->CaseType == RECONSTRUCTED_CASE)
    {
    if (this->ProcessId == 0)
      {
      ret = this->Superclass::RequestData(request, inputVector, outputVector);
      }

    this->CreateRegionArrays(output);
    this->BroadcastStructure(output, 1);

    if (this->ProcessId > 0)
      {
      output->GetFieldData()->AddArray(this->Superclass::CasePath);
      }

    this->AllReduceStatus(ret);
    this->GatherMetaData();
    this->NReadersOld = -1;
    this->MTimeOld = this->GetMTime();
    this->ShowRegionNamesOld = this->ShowRegionNames;

    this->SetErrorCode(ret ? vtkErrorCode::NoError
        : vtkErrorCode::UnknownError);
    return ret;
    }

  int nSteps = 0;
#if VTK_MAJOR_VERSION >= 6
  double requestedTimeValue(0.0);
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    requestedTimeValue
        = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    if (nSteps > 0)
      {
      outInfo->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
      }
    }
#else
  double *requestedTimeValues = NULL;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    requestedTimeValues
      = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    if (nSteps > 0)
      {
      outInfo->Set(vtkDataObject::DATA_TIME_STEPS(), requestedTimeValues, 1);
      }
    }
#endif

  const int nReaders = this->Superclass::Readers->GetNumberOfItems();
  vtkAppendCompositeDataLeaves *append = (this->CaseType == DECOMPOSED_CASE_APPENDED &&
      nReaders > 1 ? vtkAppendCompositeDataLeaves::New() : NULL);
  // append->AppendFieldDataOn();

  vtkOFFReader *reader;
  // Mark mesh as changed when the number of readers changed since
  // there is no reader to mark as mesh changed when no reader is
  // present due to some errors, e. g. no processor directory found.
  this->Superclass::MeshChanged = (nReaders != this->NReadersOld ? 1 : 0);
  this->Superclass::CurrentReaderIndex = 0;
  this->Superclass::Readers->InitTraversal();
  while ((reader
      = vtkOFFReader::SafeDownCast(this->Superclass::Readers->GetNextItemAsObject()))
      != NULL)
    {
    // even if the child readers themselves are not modified, mark
    // them as modified if "this" has been modified, since they
    // refer to the property of "this"
#if VTK_MAJOR_VERSION >= 6
    if ((nSteps > 0 && reader->SetTimeValue(requestedTimeValue))
        || this->MTimeOld != this->GetMTime())
#else
    if ((nSteps > 0 && reader->SetTimeValue(requestedTimeValues[0]))
        || this->MTimeOld != this->GetMTime())
#endif
      {
      reader->Modified();
      }
    if ((ret = reader->MakeMetaDataAtTimeStep(false)) &&
        this->CaseType == DECOMPOSED_CASE_APPENDED && nReaders > 1)
      {
      append->AddInputConnection(reader->GetOutputPort());
      }
    }

  this->GatherMetaData();

  if (this->CaseType == DECOMPOSED_CASE_APPENDED)
    {
    if (nReaders == 1)
      {
      if (ret)
        {
        reader = vtkOFFReader::SafeDownCast(
            this->Superclass::Readers->GetItemAsObject(0));
        // reader->RequestInformation() and RequestData() are called
        // for all reader instances without setting UPDATE_TIME_STEPS.
        // reader->GetExecutive()->Update() is used because
        // reader->Update() doesn't return a value.
        ret = reader->GetExecutive()->Update();
        output->ShallowCopy(reader->GetOutput());
        }
      }
    else if (nReaders > 1)
      {
      if (append->GetNumberOfInputConnections(0) != 0)
        {
        // reader->RequestInformation() and RequestData() are called
        // for all reader instances without setting UPDATE_TIME_STEPS.
        // reader->GetExecutive()->Update() is used because
        // reader->Update() doesn't return a value.
        ret = append->GetExecutive()->Update();
        output->ShallowCopy(append->GetOutput());
        }
      append->Delete();
      }

    if (!ret)
      {
      output->Initialize();
      }

    if (this->MaximumNumberOfPieces < this->NumProcesses)
      {
      this->BroadcastStructure(output, this->MaximumNumberOfPieces);
      }
    this->CreateRegionArrays(output);
    }
  else // DECOMPOSED_CASE_MULTIBLOCK
    {
    output->SetNumberOfBlocks(static_cast<unsigned int>(this->MaximumPieceId
        + 1));
    this->Superclass::CurrentReaderIndex = 0;
    this->Superclass::Readers->InitTraversal();
    for (int readerI = 0; (reader = vtkOFFReader::SafeDownCast(
        this->Superclass::Readers->GetNextItemAsObject())); readerI++)
      {
      // execute the readers
      ret *= reader->GetExecutive()->Update();
      output->SetBlock(static_cast<unsigned int>(
          this->PieceIds->GetValue(readerI)), reader->GetOutput());
      }

    if (this->NumProcesses > 1)
      {
      vtkIntArray *sendBs = vtkIntArray::New();

      // encode the multiblock structures of all pieces in this process
      for (unsigned int pieceI = 0;
          pieceI < this->PieceIds->GetNumberOfTuples(); pieceI++)
        {
        const int pieceId = this->PieceIds->GetValue(pieceI);
        sendBs->InsertNextValue(pieceId);
        vtkMultiBlockDataSet *ds = vtkMultiBlockDataSet::SafeDownCast(
            output->GetBlock(static_cast<unsigned int>(pieceId)));
        sendBs->InsertNextValue(ds->GetNumberOfBlocks());
#if VTK_MAJOR_VERSION >= 6
        vtkDataObjectTreeIterator *iter = ds->NewTreeIterator();
#else
        vtkCompositeDataIterator *iter = ds->NewIterator();
#endif
        iter->SkipEmptyNodesOff();
        iter->VisitOnlyLeavesOff();
        for (iter->InitTraversal(); !iter->IsDoneWithTraversal();
            iter->GoToNextItem())
          {
          vtkDataObject *leaf = iter->GetCurrentDataObject();
          vtkMultiBlockDataSet *mbds;
          if (leaf != NULL
              && (mbds = vtkMultiBlockDataSet::SafeDownCast(leaf)) != NULL)
            {
            sendBs->InsertNextValue(static_cast<int>(
                mbds->GetNumberOfBlocks()));
            }
          else
            {
            sendBs->InsertNextValue(-1);
            }
          }
        iter->Delete();
        }
      sendBs->Squeeze();

      // AllGather the structure. Send/receive buffers cannot be shared.
      vtkIntArray *recvBs = vtkIntArray::New();
      this->AllGatherV(sendBs, recvBs);
      sendBs->Delete();

      // decode the multiblock structure to dummy multiblock dataset
      // if the corresponding piece is not present in this process
      vtkVariantArray *mbdsArray = vtkVariantArray::New();
      vtkUnsignedIntArray *blockIArray = vtkUnsignedIntArray::New();
      vtkUnsignedIntArray *nBlocksArray = vtkUnsignedIntArray::New();
      vtkIdType bsI = 0;
      for (int pieceI = 0; pieceI < this->MaximumNumberOfPieces; pieceI++,
          bsI++)
        {
        const unsigned int pieceId
            = static_cast<unsigned int>(recvBs->GetValue(bsI++));
        vtkMultiBlockDataSet *pieceDS
            = vtkMultiBlockDataSet::SafeDownCast(output->GetBlock(pieceId));
        // if this is my block skip it; otherwise set dummy block
        if (pieceDS)
          {
          // skip the structure
          blockIArray->InsertValue(0, 0);
          int levelI = 0;
          do
            {
            unsigned int blockI = blockIArray->GetValue(levelI);
            if (blockI == 0)
              {
              nBlocksArray->InsertValue(levelI,
                  static_cast<unsigned int>(recvBs->GetValue(bsI)));
              }
            const unsigned int nBlocks = nBlocksArray->GetValue(levelI);
            for (; blockI < nBlocks; blockI++)
              {
              if (recvBs->GetValue(++bsI) >= 0)
                {
                blockIArray->SetValue(levelI++, blockI + 1);
                blockIArray->InsertValue(levelI++, 0);
                break;
                }
              }
            levelI--;
            }
          while (levelI >= 0);
          }
        else
          {
          // decode the multiblock structure and set the decoded
          // multiblock as the dummy dataset
          vtkMultiBlockDataSet *dummyBlock = output->NewInstance();
          mbdsArray->InsertValue(0, dummyBlock);
          blockIArray->InsertValue(0, 0);
          int levelI = 0;
          do
            {
            vtkMultiBlockDataSet *mbds = vtkMultiBlockDataSet::SafeDownCast(
                mbdsArray->GetValue(levelI).ToVTKObject());
            unsigned int blockI = blockIArray->GetValue(levelI);
            if (blockI == 0)
              {
              mbds->SetNumberOfBlocks(
                  static_cast<unsigned int>(recvBs->GetValue(bsI)));
              }
            const unsigned int nBlocks = mbds->GetNumberOfBlocks();
            for (; blockI < nBlocks; blockI++)
              {
              if (recvBs->GetValue(++bsI) >= 0)
                {
                vtkMultiBlockDataSet *ds = mbds->NewInstance();
                mbds->SetBlock(blockI, ds);
                ds->Delete(); // still the pointer should remain valid
                blockIArray->SetValue(levelI, blockI + 1);
                levelI++;
                mbdsArray->InsertValue(levelI, ds);
                blockIArray->InsertValue(levelI, 0);
                levelI++; // account for levelI-- at the end of the outer loop
                break;
                }
              }
            levelI--;
            }
          while (levelI >= 0);

          output->SetBlock(pieceId, dummyBlock);
          dummyBlock->Delete();
          }
        }
      mbdsArray->Delete();
      blockIArray->Delete();
      nBlocksArray->Delete();
      recvBs->Delete();
      }

    // set the processor subdirectory name as the block name. The
    // names are set to missing processor subdirectories as well.
    vtksys_ios::ostringstream os;
    for (unsigned int pieceId = 0;
        pieceId < static_cast<unsigned int>(this->MaximumPieceId + 1);
        pieceId++)
      {
      if (!vtkMultiBlockDataSet::SafeDownCast(output->GetBlock(pieceId)))
        {
        vtkMultiBlockDataSet *ds = output->NewInstance();
        output->SetBlock(pieceId, ds);
        ds->Delete();
        }
      os.str("");
      os << "processor" << pieceId;
      output->GetMetaData(pieceId)->Set(vtkCompositeDataSet::NAME(),
          os.str().c_str());
      }

    this->CreateRegionArrays(output);
    if (this->RedrawRegionNames && this->NumProcesses > 1)
      {
      if (this->ProcessId == 0)
        {
        vtkDoubleArray *recvRC = vtkDoubleArray::New();
        recvRC->SetNumberOfComponents(3);
        this->GatherV(this->RegionCentroids, recvRC);
        this->RegionCentroids->Delete();
        this->RegionCentroids = recvRC;
	vtkIntArray *recvRS = vtkIntArray::New();
	this->GatherV(this->RegionStyles, recvRS);
	this->RegionStyles->Delete();
	this->RegionStyles = recvRS;
        }
      else
        {
        this->GatherV(this->RegionCentroids, 0);
        this->GatherV(this->RegionStyles, 0);
        }
        this->GatherV(this->RegionNames);
      }
    }

  output->GetFieldData()->AddArray(this->Superclass::CasePath);

  this->Superclass::UpdateStatus();
  this->MTimeOld = this->GetMTime();
  this->ShowRegionNamesOld = this->ShowRegionNames;
  this->NReadersOld = nReaders;

  this->AllReduceStatus(ret);
  this->SetErrorCode(ret ? vtkErrorCode::NoError : vtkErrorCode::UnknownError);
  return ret;
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::CreateRegionArrays(vtkMultiBlockDataSet *output)
{
  int recreateRegionArrays = this->ShowRegionNames != this->ShowRegionNamesOld
      || (this->ShowRegionNames && this->Superclass::MeshChanged);

  if (this->NumProcesses > 1)
    {
    int recvBuffer;
    this->Controller->AllReduce(&recreateRegionArrays, &recvBuffer, 1,
        vtkCommunicator::LOGICAL_OR_OP);
    recreateRegionArrays = recvBuffer;
    }

  if (recreateRegionArrays)
    {
    this->RegionCentroids->Initialize();
    this->RegionNames->Initialize();
    this->RegionStyles->Initialize();

    if (this->ShowRegionNames)
      {
      this->MultiBlockRegionCentroids(output, 0);
      this->RegionCentroids->Squeeze();
      this->RegionNames->Squeeze();
      this->RegionStyles->Squeeze();
      }
    this->RedrawRegionNames = 1;
    }
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::MultiBlockRegionCentroids(
    vtkMultiBlockDataSet *input, const char *regionName)
{
  vtkStdString name;
  if (regionName)
    {
    name = regionName;
    }

  unsigned int nBlocks = input->GetNumberOfBlocks();
  for (unsigned int blockI = 0; blockI < nBlocks; blockI++)
    {
    vtkDataObject *dobj = input->GetBlock(blockI);
    if (!dobj)
      {
      continue;
      }

    vtkInformation *md = input->GetMetaData(blockI);

    if (dobj->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
      {
      // internal mesh
      if (md->Has(vtkCompositeDataSet::NAME()))
        {
        if (this->InternalRegionCentroid(vtkUnstructuredGrid::SafeDownCast(
            dobj)))
          {
          // if this is the single default region show its name as
          // [InternalMesh], otherwise show the region name (or
          // processor name).
          this->RegionNames->InsertNextValue(name == "" ?
              md->Get(vtkCompositeDataSet::NAME()) : regionName);
          // for the moment 0x4 stands for green color
	  this->RegionStyles->InsertNextValue(VTK_TEXT_CENTERED | 0x4);
          }
        }
      }
    else if (dobj->GetDataObjectType() == VTK_POLY_DATA)
      {
      // surface mesh
      const vtkIdType nPolys = this->AppendedPolyRegionCentroids(
          vtkPolyData::SafeDownCast(dobj));
      const vtkStdString polyName(name == "" ?
          md->Get(vtkCompositeDataSet::NAME()) : regionName);
      for (vtkIdType polyI = 0; polyI < nPolys; polyI++)
        {
        this->RegionNames->InsertNextValue(polyName);
        this->RegionStyles->InsertNextValue(VTK_TEXT_CENTERED);
        }
      }
    else if (dobj->GetDataObjectType() == VTK_MULTIBLOCK_DATA_SET)
      {
      if (md->Has(vtkCompositeDataSet::NAME()))
        {
        vtkMultiBlockDataSet *ds = vtkMultiBlockDataSet::SafeDownCast(dobj);
        vtkStdString dsName(vtkStdString(md->Get(vtkCompositeDataSet::NAME())));
        if (dsName == "[Patches]" || dsName == "[Lagrangian Particles]")
          {
          // the multiblock is a collection of patches or particles
          for (unsigned int blockJ = 0; blockJ < ds->GetNumberOfBlocks();
              blockJ++)
            {
            vtkPolyData *pd = vtkPolyData::SafeDownCast(ds->GetBlock(blockJ));
            if (pd)
              {
              const vtkIdType nPolys = this->AppendedPolyRegionCentroids(pd);
              const char *boundaryName
                  = ds->GetMetaData(blockJ)->Get(vtkCompositeDataSet::NAME());
              const vtkStdString polyName(name + (name == "" ? "" : "/")
                  + boundaryName);

              // if the boundary type is processor, set the vertical
              // text alignment depending on from and to processor
              // numbers in order to avoid overlaid thus unreadable
              // text display of the paired processor boundaries
              int regionStyle = VTK_TEXT_CENTERED;
              vtkFieldData *fd = pd->GetFieldData();
              vtkIntArray *boundaryType = vtkIntArray::SafeDownCast(
                  fd->GetArray("BoundaryType"));
              if (boundaryType && boundaryType->GetValue(0)
                  == this->Superclass::GetBoundaryTypeProcessor())
                {
                const int myProcNo = vtkIntArray::SafeDownCast(fd->GetArray(
                    "myProcNo"))->GetValue(0);
                const int neighbProcNo = vtkIntArray::SafeDownCast(fd->GetArray(
                    "neighbProcNo"))->GetValue(0);
                if (myProcNo > neighbProcNo)
                  {
                  regionStyle = VTK_TEXT_TOP;
                  }
                else
                  {
                  regionStyle = VTK_TEXT_BOTTOM;
                  }
                }
              for (vtkIdType polyI = 0; polyI < nPolys; polyI++)
                {
                this->RegionNames->InsertNextValue(polyName);
                this->RegionStyles->InsertNextValue(regionStyle);
                }
              }
            }
          }
        else if (dsName == "[Zones]")
          {
          for (unsigned int blockJ = 0; blockJ < ds->GetNumberOfBlocks();
              blockJ++)
            {
            vtkMultiBlockDataSet *dsJ = vtkMultiBlockDataSet::SafeDownCast(
                ds->GetBlock(blockJ));
            if (dsJ)
              {
              for (unsigned int blockK = 0; blockK < dsJ->GetNumberOfBlocks();
                  blockK++)
                {
                vtkDataObject *dsK = dsJ->GetBlock(blockK);
                vtkInformation *mdK = dsJ->GetMetaData(blockK);
                if (dsK && mdK->Has(vtkCompositeDataSet::NAME()))
                  {
                  if (dsK->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
                    {
                    if (this->InternalRegionCentroid(
                        vtkUnstructuredGrid::SafeDownCast(dsK)))
                      {
                      this->RegionNames->InsertNextValue(name
                          + (name == "" ? "" : "/") +
                          mdK->Get(vtkCompositeDataSet::NAME()));
                      // for the moment 0x4 stands for green color
                      this->RegionStyles->InsertNextValue(VTK_TEXT_CENTERED
                          | 0x4);
                      }
                    }
                  else if (dsK->GetDataObjectType() == VTK_POLY_DATA)
                    {
                    const vtkIdType nPolys = this->AppendedPolyRegionCentroids(
                        vtkPolyData::SafeDownCast(dsK));
                    const vtkStdString polyName(name + (name == "" ? "" : "/")
                        + mdK->Get(vtkCompositeDataSet::NAME()));
                    for (vtkIdType polyI = 0; polyI < nPolys; polyI++)
                      {
                      this->RegionNames->InsertNextValue(polyName);
                      this->RegionStyles->InsertNextValue(VTK_TEXT_CENTERED);
                      }
                    }
                  }
                }
              }
            }
          }
        else
          {
          this->MultiBlockRegionCentroids(ds,
              (name + (name == "" ? "" : "/") + dsName).c_str());
          }
        }
      }
    }
}

//-----------------------------------------------------------------------------
// Calculate centroid of given volume mesh by extracting nMaxCells
// cells and calculating weighted average of the centers with their
// volumes.
// * Averaging point cloud was not enough for meshes with varying density
// * Otherwise, could simply calculate the mid point of the bounding box?
int vtkPOFFReader::InternalRegionCentroid(vtkUnstructuredGrid *input)
{
  double centroid[3], totalVolume = 0.0;
  *centroid = centroid[1] = centroid[2] = 0.0;

  const vtkIdType nPoints = input->GetNumberOfPoints();
  const vtkIdType nCells = input->GetNumberOfCells();
  if (nPoints > 0 && nCells > 0)
    {
    // extract 100000 cells at maximum for faster calculation
    const vtkIdType nMaxCells = 100000;
    const vtkIdType interleave = (nMaxCells + nCells - 1) / nMaxCells;

    // calculate centroid weighted by the cell volume
    vtkIdList *ptIds = vtkIdList::New();
    vtkPoints *pts = vtkPoints::New();
    for (vtkIdType cellI = 0; cellI < nCells; cellI += interleave)
      {
      vtkCell *cell = input->GetCell(cellI);
      cell->Triangulate(0, ptIds, pts);
      const vtkIdType nTetras = pts->GetNumberOfPoints() / 4;
      for (vtkIdType triI = 0; triI < nTetras; triI++)
        {
        double pt0[3], pt1[3], pt2[3], pt3[3];
        const vtkIdType idx = triI * 4;
        pts->GetPoint(idx, pt0);
        pts->GetPoint(idx + 1, pt1);
        pts->GetPoint(idx + 2, pt2);
        pts->GetPoint(idx + 3, pt3);
        const double volume = vtkTetra::ComputeVolume(pt0, pt1, pt2, pt3);
        double center[3];
        vtkTetra::TetraCenter(pt0, pt1, pt2, pt3, center);
        *centroid += volume * (*center);
        centroid[1] += volume * center[1];
        centroid[2] += volume * center[2];
        totalVolume += volume;
        }
      }
    ptIds->Delete();
    pts->Delete();

    *centroid /= totalVolume;
    centroid[1] /= totalVolume;
    centroid[2] /= totalVolume;
    }

  int nRegions = 0;
  if (this->CaseType == DECOMPOSED_CASE_APPENDED && this->NumProcesses > 1)
    {
    // collect centroids and volumes calculated by each pvserver to
    // process 0 and calculate aggregated centroid
    if (this->ProcessId == 0)
      {
      double *centroids = new double [this->NumProcesses * 3];
      double *totalVolumes = new double [this->NumProcesses];
      this->Controller->Gather(centroid, centroids, 3, 0);
      this->Controller->Gather(&totalVolume, totalVolumes, 1, 0);
      totalVolume = 0.0;
      *centroid = centroid[1] = centroid[2] = 0.0;
      for (int procI = 0; procI < this->NumProcesses; procI++)
        {
        const double *centroidI = centroids + procI * 3;
        const double totalVolumeI = totalVolumes[procI];
        *centroid += totalVolumeI * (*centroidI);
        centroid[1] += totalVolumeI * centroidI[1];
        centroid[2] += totalVolumeI * centroidI[2];
        totalVolume += totalVolumeI;
        }
      delete [] centroids;
      delete [] totalVolumes;

      // if no pvserver has volume, the mesh is considered empty
      if (totalVolume > 0.0)
        {
        *centroid /= totalVolume;
        centroid[1] /= totalVolume;
        centroid[2] /= totalVolume;

        this->RegionCentroids->InsertNextTuple(centroid);
        nRegions = 1;
        }
      }
    else
      {
      this->Controller->Gather(centroid, 0, 3, 0);
      this->Controller->Gather(&totalVolume, 0, 1, 0);
      }
    }
  else
    {
    if (nPoints > 0 && nCells > 0)
      {
      this->RegionCentroids->InsertNextTuple(centroid);
      nRegions = 1;
      }
    }

  return nRegions;
}

//-----------------------------------------------------------------------------
// Compute centroid by brute-force approach when in appended
// decomposed case mode: collect all polyData to process 0; clean
// connectivity; then detect feature edges
int vtkPOFFReader::AppendedPolyRegionCentroids(vtkPolyData *input)
{
  if (this->CaseType == DECOMPOSED_CASE_APPENDED)
    {
    if (this->NumProcesses > 1)
      {
      if (this->ProcessId == 0)
        {
        // collect polyData in all pvservers and append them
        vtkPolyData **pds = new vtkPolyData * [this->NumProcesses];
        pds[0] = vtkPolyData::New();
        pds[0]->CopyStructure(input);
        for (int procI = 1; procI < this->NumProcesses; procI++)
          {
          pds[procI] = vtkPolyData::New();
          this->Controller->Receive(pds[procI], procI, 105);
          }
        vtkAppendPolyData *apd = vtkAppendPolyData::New();
        vtkPolyData *output = vtkPolyData::New();
        apd->ExecuteAppend(output, pds, this->NumProcesses);
        apd->Delete();
        vtkIdType nPolys = 0;
        for (int procI = 0; procI < this->NumProcesses; procI++)
          {
          nPolys += pds[procI]->GetNumberOfPolys();
          pds[procI]->Delete();
          }
        delete [] pds;

        int ret = 0;
        if (nPolys > 0)
          {
          // clean appended polyData
          vtkCleanPolyData *cpd = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION >= 6
          cpd->SetInputData(output);
#else
          cpd->SetInput(output);
#endif
          output->FastDelete();
          cpd->Update();
          ret = this->PolyRegionCentroids(cpd->GetOutput());
          cpd->Delete();
          }
        else
          {
          // omit cleaning if this is a point cloud
          ret = this->PolyRegionCentroids(output);
          output->Delete();
          }
        return ret;
        }
      else
        {
        // strip off field data for faster transfer
        vtkPolyData *pd = vtkPolyData::New();
        pd->CopyStructure(input);
        this->Controller->Send(pd, 0, 105);
        pd->Delete();
        return 0;
        }
      }
    else
      {
      // clean polyData appended by vtkAppendCompositeDataLeaves
      vtkCleanPolyData *cpd = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION >= 6
      cpd->SetInputData(input);
#else
      cpd->SetInput(input);
#endif
      cpd->Update();
      const int ret = this->PolyRegionCentroids(cpd->GetOutput());
      cpd->Delete();
      return ret;
      }
    }
  return this->PolyRegionCentroids(input);
}

//-----------------------------------------------------------------------------
// Split polygonal mesh by feature edges and calculate centroids with
// taking into account mainfold/non-manifold edges.  A cut-down and
// optimized version of vtkPolyDataConnectivityFilter plus
// non-manifold edge handling.
int vtkPOFFReader::PolyRegionCentroids(vtkPolyData *input)
{
  const vtkIdType numPts = input->GetNumberOfPoints();
  if (numPts < 1)
    {
    return 0; // No data to generate normals for
    }

  const vtkIdType numPolys = input->GetNumberOfPolys();
  if (numPolys < 1)
    {
    // assume lagrangian particles
    const float *points = vtkFloatArray::SafeDownCast(input->GetPoints()
        ->GetData())->GetPointer(0);
    double centroid[3];
    *centroid = centroid[1] = centroid[2] = 0.0;
    for (vtkIdType pointI = 0; pointI < numPts; pointI++)
      {
      const vtkIdType idx = 3 * pointI;
      *centroid += points[idx];
      centroid[1] += points[idx + 1];
      centroid[2] += points[idx + 2];
      }
    const double dblNPoints = static_cast<double>(numPts);
    *centroid /= dblNPoints;
    centroid[1] /= dblNPoints;
    centroid[2] /= dblNPoints;
    this->RegionCentroids->InsertNextTuple(centroid);

    return 1;
    }

  // omit feature edge detection if the polydata is a processor
  // boundary since processor boundaries tend to have lots of features
  // (but do not omit separate region detection)
  vtkIntArray *boundaryType = vtkIntArray::SafeDownCast(
      input->GetFieldData()->GetArray("BoundaryType"));
  const bool noDetectFeatures = (boundaryType && (boundaryType->GetValue(0)
      == this->Superclass::GetBoundaryTypeProcessor()
      || boundaryType->GetValue(0)
      == this->Superclass::GetBoundaryTypeInternal()));

  vtkPoints *inPts = input->GetPoints();

  // create a copy because we're modifying it
  vtkPolyData *output = vtkPolyData::New();
  output->CopyStructure(input);
  output->BuildCells(); //builds connectivity

  // create links. Build and use own point-cell links for performance.
  int *linkLoc = new int [numPts + 1];
  for (vtkIdType i = 0; i < numPts; i++)
    {
    linkLoc[i] = 0;
    }
  linkLoc[numPts] = 0;

  int *visited = new int[numPolys];
  vtkFloatArray *polyNormals = vtkFloatArray::New();
  polyNormals->SetNumberOfComponents(3);
  polyNormals->SetNumberOfTuples(numPolys);
  for (vtkIdType cellId = 0; cellId < numPolys; cellId++)
    {
    visited[cellId] = -1;
    vtkIdType npts, *pts;
    output->GetCellPoints(cellId, npts, pts);
    double n[3];
    vtkPolygon::ComputeNormal(inPts, npts, pts, n);
    polyNormals->SetTuple(cellId, n);
    // count the number of links for each cell
    for (vtkIdType pointI = 0; pointI < npts; pointI++)
      {
      linkLoc[pts[pointI]]++;
      }
    }

  // create the index by accumulating the numbers of links
  for (int pointI = 0; pointI < numPts; pointI++)
    {
    linkLoc[pointI + 1] += linkLoc[pointI];
    }

  vtkIdType *links = new vtkIdType [linkLoc[numPts]];

  // create the list
  for (vtkIdType cellI = 0; cellI < numPolys; cellI++)
    {
    vtkIdType npts, *pts;
    output->GetCellPoints(cellI, npts, pts);
    for (vtkIdType pointI = 0; pointI < npts; pointI++)
      {
      links[--linkLoc[pts[pointI]]] = cellI;
      }
    }

  const double cosAngle = cos(vtkMath::RadiansFromDegrees(60.0));

  // Traverse all cells marking those visited.  Each new search
  // starts a new connected region. Connected region grows
  // using a connected wave propagation.
  vtkIdList *wave = vtkIdList::New();
  wave->Allocate(numPts / 4 + 1, numPts);
  vtkIdList *wave2 = vtkIdList::New();
  wave2->Allocate(numPts / 4 + 1, numPts);

  vtkDoubleArray *centroids = vtkDoubleArray::New();
  centroids->SetNumberOfComponents(3);

  const float *outputPts = vtkFloatArray::SafeDownCast(output->GetPoints()
      ->GetData())->GetPointer(0);

  vtkIdType regionNumber = 0;

  vtkIdList *pointIds = vtkIdList::New();
  vtkPoints *points = vtkPoints::New();

  //visit all cells
  for (vtkIdType cellId = 0; cellId < numPolys; cellId++)
    {
    if (visited[cellId] < 0)
      {
      wave->InsertNextId(cellId);
      visited[cellId] = 1;
      vtkIdType numIds;

      double centroid[3];
      centroid[0] = centroid[1] = centroid[2] = 0.0;
      double totalArea = 0.0;
      while ((numIds = wave->GetNumberOfIds()) > 0)
        {
        for (vtkIdType i = 0; i < numIds; i++)
          {
          const vtkIdType inCellId = wave->GetId(i);

          vtkIdType npts, *pts;
          output->GetCellPoints(inCellId, npts, pts);

          double inCellCentroid[3];
          *inCellCentroid = inCellCentroid[1] = inCellCentroid[2] = 0.0;
          for (int j = 0; j < npts; j++)
            {
            const vtkIdType ptId = pts[j];
            vtkIdType idx = ptId * 3;
            *inCellCentroid += outputPts[idx];
            inCellCentroid[1] += outputPts[idx + 1];
            inCellCentroid[2] += outputPts[idx + 2];

            const vtkIdType neiPtId0 = pts[(j + 1) % npts];
            const vtkIdType neiPtId1 = pts[(j + npts - 1) % npts];
            // get the cells not where ptId but where neiPtId0 belongs
            // so that the spotted point where ptId belongs can be reused
            const vtkIdType nLinks = linkLoc[neiPtId0 + 1] - linkLoc[neiPtId0];
            const vtkIdType *linksI = links + linkLoc[neiPtId0];
            for (int linkI = 0; linkI < nLinks; linkI++)
              {
              const vtkIdType neiCellId = linksI[linkI];
              // myCell is already marked visited, so neiCellId ==
              // inCellId doesn't have to be checked
              if (visited[neiCellId] >= 0)
                {
                continue;
                }

              vtkIdType nNeiPts, *neiPts;
              output->GetCellPoints(neiCellId, nNeiPts, neiPts);
              vtkIdType spot = 0;
              for (; spot < nNeiPts; spot++)
                {
                if (neiPts[spot] == ptId)
                  {
                  // me and the neighbor cell share an edge
                  break;
                  }
                }
              if (spot == nNeiPts) // no edge is shared
                {
                continue;
                }

              if (noDetectFeatures || neiPts[(spot + 1) % nNeiPts] == neiPtId1)
                {
                // either the boundary is a processor boundary, or the
                // neighbor cell share three points, which means the
                // cell is at the reverse side of the cell. treat as
                // connected.
                wave2->InsertNextId(neiCellId);
                visited[neiCellId] = 1;
                }
              else
                {
                // this else-clause must come after the reverse side
                // check since the condition is true also to the
                // reverse side cell
                double thisNormal[3], neiNormal[3];
                polyNormals->GetTuple(inCellId, thisNormal);
                polyNormals->GetTuple(neiCellId, neiNormal);
                const double dot = vtkMath::Dot(thisNormal, neiNormal);
                if (neiPts[(spot + nNeiPts - 1) % nNeiPts] == neiPtId0)
                  {
                  // the cell is at the same side of the neighbor.
                  if (dot > cosAngle)
                    {
                    wave2->InsertNextId(neiCellId);
                    visited[neiCellId] = 1;
                    }
                  }
                else
                  {
                  // the cell is a reversed neighbor cell. mark as found
                  // only when no cell is found at the same side of the
                  // neighbor.
                  if (dot < -cosAngle)
                    {
                    wave2->InsertNextId(neiCellId);
                    visited[neiCellId] = 1;
                    }
                  }
                }//if not my reverse side
              }//for all neighbor cells of this edge in cell
            }//for all points of this cell

          vtkCell *cell = output->GetCell(inCellId);
          int cellType = cell->GetCellType();
          double area = 0.0;
          if (cellType == VTK_TRIANGLE)
            {
            area = vtkTriangle::SafeDownCast(cell)->ComputeArea();
            }
          else // VTK_QUAD || VTK_POLYGON
            {
            // vtkPolygon::ComputeArea() is unreliable so
            // triangulation is used for vtkPolygon as well
            cell->Triangulate(0, pointIds, points);
            const vtkIdType nTris = points->GetNumberOfPoints() / 3;
            for (vtkIdType triI = 0; triI < nTris; triI++)
              {
              double pt0[3], pt1[3], pt2[3];
              const vtkIdType idx = triI * 3;
              points->GetPoint(idx, pt0);
              points->GetPoint(idx + 1, pt1);
              points->GetPoint(idx + 2, pt2);
              area += vtkTriangle::TriangleArea(pt0, pt1, pt2);
              }
            }
          const double inCellWeight = area / static_cast<double>(npts);
          *centroid += *inCellCentroid * inCellWeight;
          centroid[1] += inCellCentroid[1] * inCellWeight;
          centroid[2] += inCellCentroid[2] * inCellWeight;
          totalArea += area;
          }//for all cells in this wave

        vtkIdList *tmpWave = wave;
        wave = wave2;
        wave2 = tmpWave;
        tmpWave->Reset();
        } //while wave is not empty

      regionNumber++;
      wave->Reset();
      wave2->Reset();

      *centroid /= totalArea;
      centroid[1] /= totalArea;
      centroid[2] /= totalArea;
      centroids->InsertNextTupleValue(centroid);
      }
    }

  delete [] links;
  delete [] linkLoc;
  wave->Delete();
  wave2->Delete();
  pointIds->Delete();
  points->Delete();
  delete [] visited;
  polyNormals->Delete();
  output->Delete();

  // limit the number of output regions to nMaxRegions in order to
  // avoid displaying enormous texts for staircased patches
  const vtkIdType nMaxRegions = 16;
  const vtkIdType interleave = (nMaxRegions + centroids->GetNumberOfTuples()
      - 1) / nMaxRegions;

  int nRegions = 0;
  for (vtkIdType regionI = 0; regionI < centroids->GetNumberOfTuples();
      regionI += interleave)
    {
    this->RegionCentroids->InsertNextTuple(regionI, centroids);
    nRegions++;
    }

  centroids->Delete();

  return nRegions;
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::AllReduceStatus(int &status)
{
  if (this->NumProcesses > 1)
    {
    int recvBuffer;
    this->Controller->AllReduce(&status, &recvBuffer, 1,
        vtkCommunicator::LOGICAL_AND_OP);
    status = recvBuffer;
    }
}

//-----------------------------------------------------------------------------
void vtkPOFFReader::GatherMetaData()
{
  if (this->NumProcesses > 1)
    {
    this->AllGatherV(this->Superclass::PatchDataArraySelection);
    this->AllGatherV(this->Superclass::SurfaceDataArraySelection);
    this->AllGatherV(this->Superclass::CellDataArraySelection);
    this->AllGatherV(this->Superclass::PointDataArraySelection);
    this->AllGatherV(this->Superclass::LagrangianDataArraySelection);
    // omit removing duplicated entries of LagrangianPaths as well
    // when the number of processes is 1 assuming there's no duplicate
    // entry within a process
    this->AllGatherV(this->Superclass::LagrangianPaths, true);
    }
}

//-----------------------------------------------------------------------------
// Decode and reconstruct the encoded multiblock structure. Helper
// function for BroadcastStructure().
void vtkPOFFReader::ConstructBlocks(vtkMultiBlockDataSet *ds, const int *dataTypes,
    int &leafI, vtkStringArray *names, vtkIdType &nameI)
{
  ds->SetNumberOfBlocks(dataTypes[leafI]);
  for (unsigned int blockI = 0; blockI < ds->GetNumberOfBlocks(); blockI++)
    {
    ds->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(),
        names->GetValue(nameI++).c_str());
    const int dataType = dataTypes[++leafI];
    if (dataType >= 0)
      {
      vtkDataObject *dobj = vtkDataObjectTypes::NewDataObject(dataType);
      if (dataType == VTK_MULTIBLOCK_DATA_SET)
        {
        this->ConstructBlocks(vtkMultiBlockDataSet::SafeDownCast(dobj),
            dataTypes, ++leafI, names, nameI);
        }
      ds->SetBlock(blockI, dobj);
      dobj->FastDelete();
      }
    }
}

//-----------------------------------------------------------------------------
// Broadcast the structure of a multiblock dataset. Everything has to
// be self-coded because
// vtkCommunicator::{Send,Receive}MultiBlockDataSet()s in ParaView
// 3.6.1 are broken (whereas those in 3.7-cvs have been fixed).
// The function not only reconstructs the skelton multiblocks but also
// adds empty detasets to each block so that region name labeling
// works for appended decomposed case.
void vtkPOFFReader::BroadcastStructure(vtkMultiBlockDataSet *ds, const int startProc)
{
  if (this->NumProcesses <= 1)
    {
    return;
    }

  if (this->ProcessId == 0)
    {
#if VTK_MAJOR_VERSION >= 6
    vtkDataObjectTreeIterator *iter = ds->NewTreeIterator();
#else
    vtkCompositeDataIterator *iter = ds->NewIterator();
#endif
    iter->SkipEmptyNodesOff();
    iter->VisitOnlyLeavesOff();

    // count the number of leaves and allocate buffer
    int nLeaves = 1; // in order to accommodate the number of blocks of the top
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal();
        iter->GoToNextItem())
      {
      nLeaves++;
      if (iter->GetCurrentDataObject()->GetDataObjectType()
          == VTK_MULTIBLOCK_DATA_SET)
        {
        nLeaves++;
        }
      }
    int *dataTypes = new int [nLeaves];

    // encode the multiblock structure
    vtkStringArray *names = vtkStringArray::New();
    dataTypes[0] = ds->GetNumberOfBlocks();
    int leafI = 1;
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal();
        iter->GoToNextItem(), leafI++)
      {
      names->InsertNextValue(
          iter->GetCurrentMetaData()->Has(vtkCompositeDataSet::NAME()) ?
          iter->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()) : "");
      vtkDataObject *leaf = iter->GetCurrentDataObject();
      if (leaf)
        {
        dataTypes[leafI] = leaf->GetDataObjectType();
        if (dataTypes[leafI] == VTK_MULTIBLOCK_DATA_SET)
          {
          dataTypes[++leafI] = static_cast<int>(
              vtkMultiBlockDataSet::SafeDownCast(leaf)->GetNumberOfBlocks());
          }
        }
      else
        {
        dataTypes[leafI] = -1;
        }
      }
    iter->Delete();

    // broadcast
    for (int procI = startProc; procI < this->NumProcesses; procI++)
      {
      this->Controller->Send(&nLeaves, 1, procI, 101);
      this->Controller->Send(dataTypes, nLeaves, procI, 102);
      }
    delete [] dataTypes;
    this->Broadcast(names, startProc);
    names->Delete();
    }
  else if (this->ProcessId >= startProc)
    {
    int nLeaves;
    this->Controller->Receive(&nLeaves, 1, 0, 101);
    int *dataTypes = new int [nLeaves];
    this->Controller->Receive(dataTypes, nLeaves, 0, 102);
    vtkStringArray *names = vtkStringArray::New();
    this->Broadcast(names, startProc);
    ds->Initialize();
    int leafI = 0;
    vtkIdType nameI = 0;
    this->ConstructBlocks(ds, dataTypes, leafI, names, nameI);
    delete [] dataTypes;
    names->Delete();
    }
}

//-----------------------------------------------------------------------------
// Broadcast a vtkStringArray in process 0 to all processes where
// processId is equal to or greater than startProc
void vtkPOFFReader::Broadcast(vtkStringArray *sa, int startProc)
{
  if (startProc == 0)
    {
    startProc = 1;
    }
  vtkIdType lengths[2];

  if (this->ProcessId == 0)
    {
    lengths[0] = sa->GetNumberOfTuples();
    lengths[1] = 0;
    for (int strI = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      lengths[1] += static_cast<vtkIdType>(sa->GetValue(strI).length()) + 1;
      }
    for (int procI = startProc; procI < this->NumProcesses; procI++)
      {
      this->Controller->Send(lengths, 2, procI, 103);
      }
    }
  else if (this->ProcessId >= startProc)
    {
    this->Controller->Receive(lengths, 2, 0, 103);
    }

  if (this->ProcessId == 0)
    {
    char *contents = new char [lengths[1]];
    for (int strI = 0, idx = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      const int len = static_cast<int>(sa->GetValue(strI).length()) + 1;
      memmove(contents + idx, sa->GetValue(strI).c_str(), len);
      idx += len;
      }
    for (int procI = startProc; procI < this->NumProcesses; procI++)
      {
      this->Controller->Send(contents, lengths[1], procI, 104);
      }
    delete [] contents;
    }
  else if (this->ProcessId >= startProc)
    {
    char *contents = new char [lengths[1]];
    this->Controller->Receive(contents, lengths[1], 0, 104);
    sa->Initialize();
    sa->SetNumberOfTuples(lengths[0]);
    for (int strI = 0, idx = 0; strI < lengths[0]; strI++)
      {
      sa->SetValue(strI, contents + idx);
      idx += static_cast<int>(sa->GetValue(strI).length()) + 1;
      }
    delete [] contents;
    }
}

//-----------------------------------------------------------------------------
// AllGather vtkStringArrays from and to all processes
void vtkPOFFReader::GatherV(vtkStringArray *s)
{
  vtkIdType length = 0;
  for (int strI = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    length += static_cast<vtkIdType>(s->GetValue(strI).length()) + 1;
    }
  if (this->ProcessId == 0)
    {
    vtkIdType *lengths = new vtkIdType [this->NumProcesses];
    this->Controller->Gather(&length, lengths, 1, 0);
    vtkIdType totalLength = 0;
    vtkIdType *offsets = new vtkIdType [this->NumProcesses];
    for (int procI = 0; procI < this->NumProcesses; procI++)
      {
      offsets[procI] = totalLength;
      totalLength += lengths[procI];
      }
    char *allContents = new char [totalLength], *contents = new char [length];
    for (int strI = 0, idx = 0; strI < s->GetNumberOfTuples(); strI++)
      {
      const int len = static_cast<int>(s->GetValue(strI).length()) + 1;
      memmove(contents + idx, s->GetValue(strI).c_str(), len);
      idx += len;
      }
    this->Controller->GatherV(contents, allContents, length, lengths, offsets,
        0);
    delete [] contents;
    delete [] lengths;
    delete [] offsets;
    s->Initialize();
    for (int idx = 0; idx < totalLength; idx += static_cast<int>(strlen(allContents + idx)) + 1)
      {
      s->InsertNextValue(allContents + idx);
      }
    s->Squeeze();
    delete [] allContents;
    }
  else
    {
    this->Controller->Gather(&length, 0, 1, 0);
    char *contents = new char [length];
    for (int strI = 0, idx = 0; strI < s->GetNumberOfTuples(); strI++)
      {
      const int len = static_cast<int>(s->GetValue(strI).length()) + 1;
      memmove(contents + idx, s->GetValue(strI).c_str(), len);
      idx += len;
      }
    this->Controller->GatherV(contents, 0, length, 0, 0, 0);
    delete [] contents;
    }
}

//-----------------------------------------------------------------------------
// AllGather vtkStringArrays from and to all processes
void vtkPOFFReader::AllGatherV(vtkStringArray *s, const bool unique)
{
  vtkIdType length = 0;
  for (int strI = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    length += static_cast<vtkIdType>(s->GetValue(strI).length()) + 1;
    }
  vtkIdType *lengths = new vtkIdType [this->NumProcesses];
  this->Controller->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [this->NumProcesses];
  for (int procI = 0; procI < this->NumProcesses; procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for (int strI = 0, idx = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    const int len = static_cast<int>(s->GetValue(strI).length()) + 1;
    memmove(contents + idx, s->GetValue(strI).c_str(), len);
    idx += len;
    }
  this->Controller->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  s->Initialize();
  for (int idx = 0; idx < totalLength; idx += static_cast<int>(strlen(allContents + idx)) + 1)
    {
    const char *str = allContents + idx;
    // insert only when the same string is not found
    if (s->LookupValue(str) == -1 || !unique)
      {
      s->InsertNextValue(str);
      }
    }
  s->Squeeze();
  delete [] allContents;
}

//-----------------------------------------------------------------------------
// AllGather vtkDataArraySelections from and to all processes
void vtkPOFFReader::AllGatherV(vtkDataArraySelection *s)
{
  vtkIdType length = 0;
  for (int strI = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    length += static_cast<vtkIdType>(strlen(s->GetArrayName(strI))) + 2;
    }
  vtkIdType *lengths = new vtkIdType [this->NumProcesses];
  this->Controller->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [this->NumProcesses];
  for (int procI = 0; procI < this->NumProcesses; procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for (int strI = 0, idx = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    const char *arrayName = s->GetArrayName(strI);
    contents[idx] = s->ArrayIsEnabled(arrayName);
    const int len = static_cast<int>(strlen(arrayName)) + 1;
    memmove(contents + idx + 1, arrayName, len);
    idx += len + 1;
    }
  this->Controller->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  // do not RemoveAllArray so that the previous arrays are preserved
  // s->RemoveAllArrays();
  for (int idx = 0; idx < totalLength; idx += static_cast<int>(strlen(allContents + idx + 1)) + 2)
    {
    const char *arrayName = allContents + idx + 1;
    s->AddArray(arrayName);
    if (allContents[idx] == 0)
      {
      s->DisableArray(arrayName);
      }
    else
      {
      s->EnableArray(arrayName);
      }
    }
  delete [] allContents;
}

//-----------------------------------------------------------------------------
// GatherV vtkDataArray
void vtkPOFFReader::GatherV(vtkDataArray *sendBuffer,
    vtkDataArray *recvBuffer)
{
  vtkIdType *recvLengths = new vtkIdType [this->NumProcesses];
  vtkIdType *offsets = new vtkIdType [this->NumProcesses + 1];
  const int numComponents = sendBuffer->GetNumberOfComponents();
  const vtkIdType sendLength = numComponents * sendBuffer->GetNumberOfTuples();
  this->Controller->Gather(&sendLength, recvLengths, 1, 0);
  if (this->ProcessId == 0)
    {
    offsets[0] = 0;
    for (int i = 0; i < this->NumProcesses; i++)
      {
      offsets[i+1] = offsets[i] + recvLengths[i];
      }
    recvBuffer->SetNumberOfComponents(numComponents);
    recvBuffer->SetNumberOfTuples(offsets[this->NumProcesses] / numComponents);
    }
  this->Controller->GetCommunicator()->GatherVVoidArray(
      sendBuffer->GetVoidPointer(0), (recvBuffer
      ? recvBuffer->GetVoidPointer(0) : 0), sendLength, recvLengths, offsets,
      sendBuffer->GetDataType(), 0);
  delete [] offsets;
  delete [] recvLengths;
}

//-----------------------------------------------------------------------------
// AllGatherV vtkDataArray
// the function is required because
// vtkMultiProcessController::AllGatherV(vtkDataArray *, vtkDataArray *)
// is not in PV 3.8 and VTK 5.6 but is only in PV 3.9-git and VTK 5.7-git.
void vtkPOFFReader::AllGatherV(vtkDataArray *sendBuffer,
    vtkDataArray *recvBuffer)
{
  vtkIdType *recvLengths = new vtkIdType [this->NumProcesses];
  vtkIdType *offsets = new vtkIdType [this->NumProcesses + 1];
  const int numComponents = sendBuffer->GetNumberOfComponents();
  const vtkIdType sendLength = numComponents * sendBuffer->GetNumberOfTuples();
  this->Controller->AllGather(&sendLength, recvLengths, 1);
  offsets[0] = 0;
  for (int i = 0; i < this->NumProcesses; i++)
    {
    offsets[i + 1] = offsets[i] + recvLengths[i];
    }
  recvBuffer->SetNumberOfComponents(numComponents);
  recvBuffer->SetNumberOfTuples(offsets[this->NumProcesses] / numComponents);
  this->Controller->GetCommunicator()->AllGatherVVoidArray(
      sendBuffer->GetVoidPointer(0), recvBuffer->GetVoidPointer(0),
      sendLength, recvLengths, offsets, sendBuffer->GetDataType());
  delete [] offsets;
  delete [] recvLengths;
}
