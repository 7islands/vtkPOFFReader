/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// Thanks to Terry Jordan of SAIC at the National Energy
// Technology Laboratory who developed this class.
// Please address all comments to Terry Jordan (terry.jordan@sa.netl.doe.gov)
//
// Token-based FoamFile format lexer / parser, performance / stability
// / compatibility enhancements, gzipped file support, lagrangian
// field support, variable timestep support, builtin cell-to-point
// filter, pointField support, polyhedron decomposition support, OF
// 1.5 extended format support, multiregion support, old mesh format
// support, parallelization support for decomposed cases in
// conjunction with vtkPOFFReader, et. al. by Takuya Oshima of
// Niigata University, Japan (oshima@eng.niigata-u.ac.jp)
//
// * GUI Based selection of mesh regions and fields available in the case
// * Minor bug fixes / Strict memory allocation checks
// * Minor performance enhancements
// by Philippose Rajan (sarith@rocketmail.com)
//
// Disclaimer:
// OPENFOAM(R) is a registered trade mark of OpenCFD Limited, the
// producer of the OpenFOAM software and owner of the OPENFOAM(R) and
// OpenCFD(R) trade marks. This code is not approved or endorsed by
// OpenCFD Limited.


// Hijack the CRC routine of zlib to omit CRC check for gzipped files
// (on OSes other than Windows where the mechanism doesn't work due
// to pre-bound DLL symbols) if set to 1, or not (set to 0). Affects
// performance by about 3% - 4%.
#define VTK_FOAMFILE_OMIT_CRCCHECK 0

// The input/output buffer sizes for zlib in bytes.
#define VTK_FOAMFILE_INBUFSIZE (16384)
#define VTK_FOAMFILE_OUTBUFSIZE (131072)
#define VTK_FOAMFILE_INCLUDE_STACK_SIZE (10)

// Avoid a locale problem where period is not interpreted as decimal
// point when ParaView is compiled with Qt 4.5 and run under certain
// locales (e.g. de_DE, fr_FR)
#ifndef VTK_FOAMFILE_LOCALE_WORKAROUND
#define VTK_FOAMFILE_LOCALE_WORKAROUND 1
#endif

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#if VTK_FOAMFILE_OMIT_CRCCHECK
#define ZLIB_INTERNAL
#endif

// for possible future extension of linehead-aware directives
#define VTK_FOAMFILE_RECOGNIZE_LINEHEAD 0

#include "vtkOFFDevReader.h"

#if !defined(VTK_FOAMFILE_HAVE_REGEX)
#include "vtksys/RegularExpression.hxx"
#endif
#include <vtkstd/vector>
#include "vtksys/DateStamp.h"
#include "vtksys/SystemTools.hxx"
#include <vtksys/ios/sstream>
#include "vtk_zlib.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCharArray.h"
#include "vtkCollection.h"
#include "vtkConvexPointSet.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFloatArray.h"
#include "vtkHexahedron.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkPyramid.h"
#include "vtkQuad.h"
#include "vtkSortDataArray.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVertex.h"
#include "vtkWedge.h"

#if !(defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__))
// for getpwnam() / getpwuid()
#include <sys/types.h>
#include <pwd.h>
// for regcomp() / regexec() / regfree() / regerror()
#include <regex.h>
// for getuid()
#include <unistd.h>
#endif
// for fabs()
#include <math.h>
// for isalnum() / isspace() / isdigit()
#include <ctype.h>

// If 1, Use VTK_POLYHEDRON in order to represent polyhedral cells
// without decomposition. If 0, use VTK_CONVEX_POINT_SET instead
#if vtksys_DATE_STAMP_FULL >= 20100525
#define VTK_FOAMFILE_USE_VTK_POLYHEDRON 1
#else
#define VTK_FOAMFILE_USE_VTK_POLYHEDRON 0
#endif

// If 1, reorder symmTensor components so that they match the
// component label display in ParaView. Otherwise relabel the
// component panel display in ParaView so that it matches the OpenFOAM
// symmTensor order.
#if vtksys_DATE_STAMP_FULL < 20100429
#define VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS 1
#else
#define VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS 0
#endif

// avoid name crashes with the builtin reader
#define vtkFoamArrayVector vtkNewFoamArrayVector
#define vtkFoamError vtkNewFoamError
#define vtkFoamToken vtkNewFoamToken
#define vtkFoamFileStack vtkNewFoamFileStack
#define vtkFoamFile vtkNewFoamFile
#define vtkFoamIOobject vtkNewFoamIOobject
#define vtkFoamReadValue vtkNewFoamReadValue
#define vtkFoamEntryValue vtkNewFoamEntryValue
#define vtkFoamEntry vtkNewFoamEntry
#define vtkFoamDict vtkNewFoamDict

#if VTK_FOAMFILE_OMIT_CRCCHECK
uLong ZEXPORT crc32(uLong, const Bytef *, uInt)
{ return 0; }
#endif

vtkStandardNewMacro(vtkOFFReader);

// forward declarations
template <typename T> struct vtkFoamArrayVector
  : public vtkstd::vector<T *>
{
private:
  typedef vtkstd::vector<T *> Superclass;

public:
  ~vtkFoamArrayVector()
  {
    for(size_t arrayI = 0; arrayI < Superclass::size(); arrayI++)
      {
      if(Superclass::operator[](arrayI))
        {
        Superclass::operator[](arrayI)->Delete();
        }
      }
  }
};

typedef vtkFoamArrayVector<vtkIdList> vtkFoamIdListVector;
typedef vtkFoamArrayVector<vtkIntArray> vtkFoamIntArrayVector;
typedef vtkFoamArrayVector<vtkFloatArray> vtkFoamFloatArrayVector;
struct vtkFoamIntVectorVector;

struct vtkFoamError;
struct vtkFoamToken;
struct vtkFoamFileStack;
struct vtkFoamFile;
struct vtkFoamIOobject;
template <typename T> struct vtkFoamReadValue;
struct vtkFoamEntryValue;
struct vtkFoamEntry;
struct vtkFoamDict;

//-----------------------------------------------------------------------------
// class vtkOFFReaderPrivate
// the reader core of vtkOFFReader
class vtkOFFReaderPrivate : public vtkObject
{
public:
  static vtkOFFReaderPrivate *New();
  vtkTypeMacro(vtkOFFReaderPrivate, vtkObject);
  void PrintSelf(ostream &, vtkIndent);

  vtkDoubleArray *GetTimeValues()
    {return this->TimeValues;}
  vtkStringArray *GetTimeNames()
    {return this->TimeNames;}
  vtkGetMacro(TimeStep, int);
  vtkSetMacro(TimeStep, int);
  const vtkStdString &GetRegionName() const
    {return this->RegionName;}

  // gather timestep information
  bool MakeInformationVector(const vtkStdString &, const vtkStdString &,
      const vtkStdString &, vtkOFFReader *);
  // read mesh/fields and create dataset
  int RequestData(vtkMultiBlockDataSet *, bool, bool, bool, bool);
  void SetTimeValue(const double);
  int MakeMetaDataAtTimeStep(vtkStringArray *, vtkStringArray *,
      vtkStringArray *, vtkStringArray *, const bool);
  void SetupInformation(const vtkStdString &, const vtkStdString &,
      const vtkStdString &, vtkOFFReaderPrivate *);
  static int GetBoundaryTypeProcessor();
  static int GetBoundaryTypeInternal();

private:
  struct vtkFoamBoundaryEntry
    {
    enum bt
      {
      // each value must be in 2^n for bitwise operations
      GEOMETRICAL = 0, // symmetryPlane, wedge, cyclic, empty, etc.
      PHYSICAL = 1, // patch, wall
      PROCESSOR = 2, // processor
      INTERNAL = 4 // internal faces
      };
    vtkStdString BoundaryName;
    int NFaces, StartFace, AllBoundariesStartFace, MyProcNo, NeighbProcNo;
    bool IsActive;
    bt BoundaryType;
    };

  struct vtkFoamBoundaryDict : public vtkstd::vector<vtkFoamBoundaryEntry>
    {
    // we need to keep the path to time directory where the current mesh
    // is read from, since boundaryDict may be accessed multiple times
    // at a timestep for patch selections
    vtkStdString TimeDir;

    // retain face number ranges where myProcNo is larger than
    // neighbProcNo so that polyhedral cell faces on a processor
    // boundary can be decomposed properly
    vtkstd::vector<int> UpperProcRanges;
    int LowestUpperProcFaceNo;

    // processor boundaries used to reconstruct processor faces when
    // in appended decomposed case mode.
    vtkstd::vector<int> ProcBoundaries;
    };

  vtkOFFReader *Parent;

  // case and region
  vtkStdString CasePath;
  vtkStdString RegionName;
  vtkStdString ProcessorName;

  // time information
  vtkDoubleArray *TimeValues;
  int TimeStep;
  int TimeStepOld;
  vtkStringArray *TimeNames;

  int InternalMeshSelectionStatus;
  int InternalMeshSelectionStatusOld;

  int SurfaceMeshSelectionStatus;
  int SurfaceMeshSelectionStatusOld;

  // filenames / directories
  vtkStringArray *VolFieldFiles;
  vtkStringArray *SurfaceFieldFiles;
  vtkStringArray *PointFieldFiles;
  vtkStringArray *LagrangianFieldFiles;
  vtkStringArray *PolyMeshPointsDir;
  vtkStringArray *PolyMeshFacesDir;

  // for mesh construction
  vtkIdType NumCells;
  vtkIdType NumPoints;
  vtkIdType NumInternalFaces;
  vtkIntArray *FaceOwner;
  vtkFoamIntVectorVector *ProcessorFaces;

  // for cell-to-point interpolation
  vtkPolyData *AllBoundaries;
  vtkIntArray *AllBoundariesPointMap;
  vtkIntArray *InternalPoints;

  // for caching mesh
  vtkUnstructuredGrid *InternalMesh;
  vtkPolyData *SurfaceMesh;
  vtkMultiBlockDataSet *BoundaryMesh;
  vtkFoamIntArrayVector *BoundaryPointMap; // use IntArray for LookupValue()
  vtkFoamBoundaryDict BoundaryDict;
  vtkMultiBlockDataSet *LagrangianMesh;
  vtkMultiBlockDataSet *PointZoneMesh;
  vtkMultiBlockDataSet *FaceZoneMesh;
  vtkMultiBlockDataSet *CellZoneMesh;
#if 0
  vtkFoamFloatArrayVector *ReciprocalDelta;
#endif

  // for polyhedra handling
  int NumTotalAdditionalCells;
  vtkIntArray *AdditionalCellIds;
  vtkIntArray *NumAdditionalCells;
  vtkFoamIdListVector *AdditionalCellPoints;

  // region names displayed in the selection list
  static const char *InternalMeshIdentifier;
  static const char *SurfaceMeshIdentifier;
  static const char *ProcessorPatchesIdentifier;

  // constructor and destructor are kept private
  vtkOFFReaderPrivate();
  ~vtkOFFReaderPrivate();

  // not implemented.
  vtkOFFReaderPrivate(const vtkOFFReaderPrivate &);
  void operator=(const vtkOFFReaderPrivate &);

  // clear mesh construction
  void ClearInternalMeshes();
  void ClearBoundaryMeshes();
  void ClearLagrangianMeshes();
  void ClearMeshes();

  vtkStdString RegionPath() const
    {return (this->RegionName == "" ? "" : "/") + this->RegionName;}
  vtkStdString TimePath(const int timeI) const
    {return this->CasePath + this->TimeNames->GetValue(timeI);}
  vtkStdString TimeRegionPath(const int timeI) const
    {return this->TimePath(timeI) + this->RegionPath();}
  vtkStdString TimeRegionMeshPath(vtkStringArray *dir, const int timeI) const
    {return this->CasePath + dir->GetValue(timeI) + this->RegionPath()
    + "/polyMesh/";}
  vtkStdString CurrentTimePath() const
    {return this->TimePath(this->TimeStep);}
  vtkStdString CurrentTimeRegionPath() const
    {return this->TimeRegionPath(this->TimeStep);}
  vtkStdString CurrentTimeRegionMeshPath(vtkStringArray *dir) const
    {return this->TimeRegionMeshPath(dir, this->TimeStep);}
  vtkStdString RegionPrefix() const
    {return this->RegionName + (this->RegionName == "" ? "" : "/");}

  // search time directories for mesh
  void AppendMeshDirToArray(vtkStringArray *, const vtkStdString &,
      const vtkStdString &, const int);
  void PopulatePolyMeshDirArrays();

  // search a time directory for field objects
  void GetFieldNames(const vtkStdString &, const bool, vtkStringArray *,
      vtkStringArray *, vtkStringArray *);
  void SortFieldFiles(vtkStringArray *, vtkStringArray *, vtkStringArray *);
  void LocateLagrangianClouds(vtkStringArray *, const vtkStdString &);

  // read controlDict
  bool ListTimeDirectoriesByControlDict(vtkFoamDict *dict);
  bool ListTimeDirectoriesByInstances();

  // read mesh files
  vtkFloatArray* ReadPointsFile();
  vtkFoamIntVectorVector* ReadFacesFile (const vtkStdString &);
  vtkFoamIntVectorVector* ReadOwnerNeighborFiles(const vtkStdString &,
      vtkFoamIntVectorVector *);
  bool CheckFacePoints(vtkFoamIntVectorVector *);

  // create mesh
  void InsertCellsToGrid(vtkUnstructuredGrid *, const vtkFoamIntVectorVector *,
      const vtkFoamIntVectorVector *, vtkFloatArray *, vtkIdTypeArray *,
      vtkIntArray *);
  vtkUnstructuredGrid *MakeInternalMesh(const vtkFoamIntVectorVector *,
      const vtkFoamIntVectorVector *, vtkFloatArray *);
  vtkPolyData *MakeSurfaceMesh(const vtkStdString &,
      const vtkFoamIntVectorVector *, vtkFloatArray *);

  bool InsertFacesToGrid(vtkPolyData *, const vtkFoamIntVectorVector *, int,
      int, vtkIntArray *, vtkIdList *, const int *, const bool);
  template <typename T1, typename T2> bool ExtendArray(T1 *, const int);
  vtkMultiBlockDataSet* MakeBoundaryMesh(const vtkFoamIntVectorVector *,
      vtkFloatArray *);
  void SetBlockName(vtkMultiBlockDataSet *, unsigned int, const char *);
  void TruncateFaceOwner();
#if 0
  void CalculateReciprocalDelta(const vtkFoamIntVectorVector *);
#endif

  // move additional points for decomposed cells
  vtkPoints *MoveInternalMesh(vtkUnstructuredGrid *, vtkFloatArray *);
  int MoveSurfaceMesh(vtkPolyData *internalMesh, vtkFloatArray *pointArray);
  void MoveBoundaryMesh(vtkMultiBlockDataSet *, vtkFloatArray *);

  // cell-to-point interpolator
  void InterpolateCellToPoint(vtkFloatArray *, vtkFloatArray *, vtkPointSet *,
      vtkIntArray *, const int);

  // read and create cell/point fields
  void ConstructDimensions(vtkStdString *, vtkFoamDict *);
  bool ReadFieldFile(vtkFoamIOobject *, vtkFoamDict *, const vtkStdString &,
      vtkDataArraySelection *);
  vtkFloatArray *FillField(vtkFoamEntry *, vtkFoamEntry *, int, vtkFoamIOobject *,
      const vtkStdString &, const bool);
  void GetVolFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
      const vtkStdString &);
  void GetSurfaceFieldAtTimeStep(vtkPolyData *, vtkMultiBlockDataSet *,
      const vtkStdString &);

  void GetPointFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
      const vtkStdString &);
  void AddArrayToFieldData(vtkDataSetAttributes *, vtkDataArray *,
      const vtkStdString &);

  // create lagrangian mesh/fields
  vtkMultiBlockDataSet *MakeLagrangianMesh();

  // create point/face/cell zones
  vtkFoamDict *GatherBlocks(const char *, const int);
  bool GetPointZoneMesh(vtkMultiBlockDataSet *, vtkPoints *);
  bool GetFaceZoneMesh(vtkMultiBlockDataSet *, const vtkFoamIntVectorVector *,
      vtkPoints *);
  bool GetCellZoneMesh(vtkMultiBlockDataSet *, const vtkFoamIntVectorVector *,
      const vtkFoamIntVectorVector *, vtkPoints *);
};

vtkStandardNewMacro(vtkOFFReaderPrivate);

//-----------------------------------------------------------------------------
// struct vtkFoamIntVectorVector
struct vtkFoamIntVectorVector
{
private:
  vtkIntArray *Indices, *Body;

public:
  ~vtkFoamIntVectorVector()
  {
    this->Indices->Delete();
    this->Body->Delete();
  }

  vtkFoamIntVectorVector(const vtkFoamIntVectorVector &ivv) :
    Indices(ivv.Indices), Body(ivv.Body)
  {
    this->Indices->Register(0); // vtkDataArrays do not have ShallowCopy
    this->Body->Register(0);
  }
  vtkFoamIntVectorVector() :
    Indices(vtkIntArray::New()), Body(vtkIntArray::New())
  {
  }
  vtkFoamIntVectorVector(const int nElements, const int bodyLength) :
    Indices(vtkIntArray::New()), Body(vtkIntArray::New())
  {
    this->Indices->SetNumberOfValues(nElements + 1);
    this->Body->SetNumberOfValues(bodyLength);
  }

  // note that vtkIntArray::Resize() allocates (current size + new
  // size) bytes if current size < new size until 2010-06-27
  // cf. commit c869c3d5875f503e757b64f2fd1ec349aee859bf
  void ResizeBody(const int bodyLength)
  {
    this->Body->Resize(bodyLength);
  }

  int *WritePointer(const int i, const int bodyI, const int number)
  {
    return this->Body->WritePointer(*this->Indices->GetPointer(i) = bodyI,
        number);
  }

  int *SetIndex(const int i, const int bodyI)
  {
    return this->Body->GetPointer(*this->Indices->GetPointer(i) = bodyI);
  }
  void SetValue(const int bodyI, int value)
  {
    this->Body->SetValue(bodyI, value);
  }
  void InsertValue(const int bodyI, int value)
  {
    this->Body->InsertValue(bodyI, value);
  }
  const int *operator[](const int i) const
  {
    return this->Body->GetPointer(this->Indices->GetValue(i));
  }
  int GetSize(const int i) const
  {
    return this->Indices->GetValue(i + 1) - this->Indices->GetValue(i);
  }
  int GetNumberOfElements() const
  {
    return this->Indices->GetNumberOfTuples() - 1;
  }
  vtkIntArray *GetIndices()
  {
    return this->Indices;
  }
  vtkIntArray *GetBody()
  {
    return this->Body;
  }
};

//-----------------------------------------------------------------------------
// class vtkFoamError
// class for exception-carrying object
struct vtkFoamError : public vtkStdString
{
private:
  typedef vtkStdString Superclass;

public:
  vtkFoamError() :
    vtkStdString()
  {
  }
  vtkFoamError(const vtkFoamError& e) :
    vtkStdString(e)
  {
  }
  ~vtkFoamError()
  {
  }
  // a super-easy way to make use of operator<<()'s defined in
  // vtksys_ios::ostringstream class
  template <class T> vtkFoamError& operator<<(const T& t)
  {
    vtksys_ios::ostringstream os;
    os << t;
    this->Superclass::operator+=(os.str());
    return *this;
  }
};

//-----------------------------------------------------------------------------
// class vtkFoamToken
// token class which also works as container for list types
// - a word token is treated as a string token for simplicity
// - handles only atomic types. Handling of list types are left to the
//   derived classes.
struct vtkFoamToken
{
public:
  enum tokenType
    {
    // undefined type
    UNDEFINED,
    // atomic types
    PUNCTUATION, LABEL, SCALAR, WORD, STRING, IDENTIFIER,
    // vtkObject-derived list types
    STRINGLIST, LABELLIST, SCALARLIST, VECTORLIST,
    // original list types
    LABELLISTLIST, ENTRYVALUELIST, EMPTYLIST, DICTIONARY,
    // error state
    TOKEN_ERROR
    };

protected:
  tokenType Type;
  union
  {
    char Char;
    int Int;
    double Double;
    vtkStdString* String;
    vtkObjectBase *VtkObjectPtr;
    // vtkObject-derived list types
    vtkIntArray *LabelListPtr;
    vtkFloatArray *ScalarListPtr, *VectorListPtr;
    vtkStringArray *StringListPtr;
    // original list types
    vtkFoamIntVectorVector *LabelListListPtr;
    vtkstd::vector<vtkFoamEntryValue*> *EntryValuePtrs;
    vtkFoamDict *DictPtr;
  };

  void Clear()
  {
    if (this->Type == WORD || this->Type == STRING || this->Type == IDENTIFIER)
      {
      delete this->String;
      }
  }

  void AssignData(const vtkFoamToken& value)
  {
    switch (value.Type)
      {
      case PUNCTUATION:
        this->Char = value.Char;
        break;
      case LABEL:
        this->Int = value.Int;
        break;
      case SCALAR:
        this->Double = value.Double;
        break;
      case WORD:
      case STRING:
      case IDENTIFIER:
        this->String = new vtkStdString(*value.String);
        break;
        // required to suppress the 'enumeration value not handled' warning by
        // g++ when compiled with -Wall
      default:
        break;
      }
  }

public:
  vtkFoamToken() :
    Type(UNDEFINED)
  {
  }
  vtkFoamToken(const vtkFoamToken& value) :
    Type(value.Type)
  {
    this->AssignData(value);
  }
  virtual ~vtkFoamToken()
  {
    this->Clear();
  }

  tokenType GetType() const
  {
    return this->Type;
  }

  template <typename T> bool Is() const;
  template <typename T> T To() const;
#if defined(_MSC_VER)
  // workaround for Win32-64ids-nmake70
  VTK_TEMPLATE_SPECIALIZE bool Is<int>() const;
  VTK_TEMPLATE_SPECIALIZE bool Is<float>() const;
  VTK_TEMPLATE_SPECIALIZE bool Is<double>() const;
  VTK_TEMPLATE_SPECIALIZE int To<int>() const;
  VTK_TEMPLATE_SPECIALIZE float To<float>() const;
  VTK_TEMPLATE_SPECIALIZE double To<double>() const;
#endif

  // workaround for SunOS-CC5.6-dbg
  int ToInt() const
  {
    return this->Int;
  }

  // workaround for SunOS-CC5.6-dbg
  float ToFloat() const
  {
    return this->Type == LABEL ? this->Int : this->Double;
  }

  bool IsWordOrString() const
  {
    return this->Type == WORD || this->Type == STRING;
  }

  const vtkStdString &ToStdString() const
  {
    return *this->String;
  }

  void SetBad()
  {
    this->Clear();
    this->Type = TOKEN_ERROR;
  }
  void SetWord(const vtkStdString &wordString)
  {
    this->operator=(wordString);
    this->Type = WORD;
  }
  void SetIdentifier(const vtkStdString& idString)
  {
    this->operator=(idString);
    this->Type = IDENTIFIER;
  }

  void operator=(const char value)
  {
    this->Clear();
    this->Type = PUNCTUATION;
    this->Char = value;
  }
  void operator=(const int value)
  {
    this->Clear();
    this->Type = LABEL;
    this->Int = value;
  }
  void operator=(const double value)
  {
    this->Clear();
    this->Type = SCALAR;
    this->Double = value;
  }
  void operator=(const char *value)
  {
    this->Clear();
    this->Type = STRING;
    this->String = new vtkStdString(value);
  }
  void operator=(const vtkStdString& value)
  {
    this->Clear();
    this->Type = STRING;
    this->String = new vtkStdString(value);
  }
  void operator=(const vtkFoamToken& value)
  {
    this->Clear();
    this->Type = value.Type;
    this->AssignData(value);
  }
  bool operator==(const char value) const
  {
    return this->Type == PUNCTUATION && this->Char == value;
  }
  bool operator==(const int value) const
  {
    return this->Type == LABEL && this->Int == value;
  }
  bool operator==(const vtkStdString& value) const
  {
    return this->Type == WORD && *this->String == value;
  }
  bool operator!=(const char value) const
  {
    return !this->operator==(value);
  }
  bool operator!=(const vtkStdString& value) const
  {
    return !this->operator==(value);
  }

  friend vtksys_ios::ostringstream& operator<<(vtksys_ios::ostringstream& str,
      const vtkFoamToken& value)
  {
    switch (value.GetType())
      {
      case TOKEN_ERROR:
        str << "badToken (an unexpected EOF?)";
        break;
      case PUNCTUATION:
        str << value.Char;
        break;
      case LABEL:
        str << value.Int;
        break;
      case SCALAR:
        str << value.Double;
        break;
      case WORD:
      case STRING:
        str << *value.String;
        break;
      case IDENTIFIER:
        str << "$" << *value.String;
        break;
      // required to suppress the 'enumeration value not handled' warning by
      // g++ when compiled with -Wall
      default:
        break;
      }
    return str;
  }
};

VTK_TEMPLATE_SPECIALIZE inline bool vtkFoamToken::Is<int>() const
{
  return this->Type == LABEL;
}

VTK_TEMPLATE_SPECIALIZE inline bool vtkFoamToken::Is<float>() const
{
  return this->Type == LABEL || this->Type == SCALAR;
}

VTK_TEMPLATE_SPECIALIZE inline bool vtkFoamToken::Is<double>() const
{
  return this->Type == SCALAR;
}

VTK_TEMPLATE_SPECIALIZE inline int vtkFoamToken::To<int>() const
{
  return this->Int;
}

VTK_TEMPLATE_SPECIALIZE inline float vtkFoamToken::To<float>() const
{
  return this->Type == LABEL ? this->Int : this->Double;
}

VTK_TEMPLATE_SPECIALIZE inline double vtkFoamToken::To<double>() const
{
  return this->Type == LABEL ? this->Int : this->Double;
}

//-----------------------------------------------------------------------------
// class vtkFoamFileStack
// list of variables that have to be saved when a file is included.
struct vtkFoamFileStack
{
protected:
  vtkStdString FileName;
  FILE *File;
  bool IsCompressed;
  z_stream Z;
  int ZStatus;
  int LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
  bool WasNewline;
#endif

  // buffer pointers. using raw pointers for performance reason.
  unsigned char *Inbuf;
  unsigned char *Outbuf;
  unsigned char *BufPtr;
  unsigned char *BufEndPtr;

  vtkFoamFileStack() :
    FileName(), File(NULL), IsCompressed(false), ZStatus(Z_OK), LineNumber(0),
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        WasNewline(true),
#endif
        Inbuf(NULL), Outbuf(NULL), BufPtr(NULL), BufEndPtr(NULL)
  {
    this->Z.zalloc = Z_NULL;
    this->Z.zfree = Z_NULL;
    this->Z.opaque = Z_NULL;
  }

  void Reset()
  {
    // this->FileName = "";
    this->File = NULL;
    this->IsCompressed = false;
    // this->ZStatus = Z_OK;
    this->Z.zalloc = Z_NULL;
    this->Z.zfree = Z_NULL;
    this->Z.opaque = Z_NULL;
    // this->LineNumber = 0;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    this->WasNewline = true;
#endif

    this->Inbuf = NULL;
    this->Outbuf = NULL;
    // this->BufPtr = NULL;
    // this->BufEndPtr = NULL;
  }

public:
  virtual ~vtkFoamFileStack()
  {
  }

  const vtkStdString& GetFileName() const
  {
    return this->FileName;
  }
  int GetLineNumber() const
  {
    return this->LineNumber;
  }
};

//-----------------------------------------------------------------------------
// class vtkFoamFile
// read and tokenize the input.
struct vtkFoamFile : public vtkFoamFileStack
{
private:
  typedef vtkFoamFileStack Superclass;

public:
  // #inputMode values
  enum inputModes
  {
    INPUT_MODE_PROTECT, INPUT_MODE_MERGE, INPUT_MODE_OVERWRITE, INPUT_MODE_WARN,
    INPUT_MODE_ERROR
  };

private:
  inputModes InputMode;

  // inclusion handling
  vtkFoamFileStack *Stack[VTK_FOAMFILE_INCLUDE_STACK_SIZE];
  int StackI;
  vtkStdString CasePath;

  // declare and define as private
  vtkFoamFile();
  bool InflateNext(unsigned char *buf, int requestSize);
  int NextTokenHead();
  // hacks to keep exception throwing / recursive codes out-of-line to make
  // putBack(), getc() and readExpecting() inline expandable
  void ThrowDuplicatedPutBackException();
  void ThrowUnexpectedEOFException();
  void ThrowUnexpectedNonNumberTokenException(const vtkFoamToken &token);
  void ThrowUnexpectedTokenException(const char, const int c);
  int ReadNext();

  void PutBack(const int c)
  {
    if (--this->Superclass::BufPtr < this->Superclass::Outbuf)
      {
      this->ThrowDuplicatedPutBackException();
      }
    *this->Superclass::BufPtr = c;
  }

  // get a character
  int Getc()
  {
    return this->Superclass::BufPtr == this->Superclass::BufEndPtr ? this->ReadNext()
        : *this->Superclass::BufPtr++;
  }

  vtkFoamError StackString()
  {
    vtksys_ios::ostringstream os;
    if (this->StackI > 0)
      {
      os << "\n included";

      for (int stackI = this->StackI - 1; stackI >= 0; stackI--)
        {
        os << " from line " << this->Stack[stackI]->GetLineNumber() << " of "
            << this->Stack[stackI]->GetFileName() << "\n";
        }
      os << ": ";
      }
    return vtkFoamError() << os.str();
  }

  bool CloseIncludedFile()
  {
    if (this->StackI == 0)
      {
      return false;
      }
    this->Clear();
    this->StackI--;
    // use the default bitwise assignment operator
    this->Superclass::operator=(*this->Stack[this->StackI]);
    delete this->Stack[this->StackI];
    return true;
  }

  void Clear()
  {
    if (this->Superclass::IsCompressed)
      {
      inflateEnd(&this->Superclass::Z);
      }

    delete [] this->Superclass::Inbuf;
    delete [] this->Superclass::Outbuf;
    this->Superclass::Inbuf = this->Superclass::Outbuf = NULL;

    if (this->Superclass::File)
      {
      fclose(this->Superclass::File);
      this->Superclass::File = NULL;
      }
    // don't reset the line number so that the last line number is
    // retained after close
    // lineNumber_ = 0;
  }

  const vtkStdString ExtractPath(const vtkStdString& path) const
  {
#if defined(_WIN32)
    const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
    const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
    const vtkStdString::size_type pos = path.find_last_of(pathFindSeparator);
    return pos == vtkStdString::npos ? vtkStdString(".") + pathSeparator
        : path.substr(0, pos + 1);
  }

  const vtkStdString ExtractCaseName() const
  {
#if defined(_WIN32)
    const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
    const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
    const vtkStdString::size_type len = this->CasePath.length();
    if (len == 0 || (len == 1
      && vtkStdString::npos != this->CasePath.find_first_of(pathFindSeparator)))
      {
      return vtkStdString();
      }
    const vtkStdString::size_type pos = this->CasePath.find_last_of(pathFindSeparator, len - 2);
    return pos == vtkStdString::npos ? this->CasePath.substr(0, len - 1)
      : this->CasePath.substr(pos + 1, len - 2 - pos);
  }

public:
  vtkFoamFile(const vtkStdString& casePath) :
    vtkFoamFileStack(), InputMode(INPUT_MODE_MERGE), StackI(0),
    CasePath(casePath)
  {
  }
  ~vtkFoamFile()
  {
    this->Close();
  }

  inputModes GetInputMode() const
  {
    return this->InputMode;
  }
  const vtkStdString GetCasePath() const
  {
    return this->CasePath;
  }
  const vtkStdString GetFilePath() const
  {
    return this->ExtractPath(this->FileName);
  }

  vtkStdString ExpandPath(const vtkStdString& pathIn,
      const vtkStdString& defaultPath)
  {
    vtkStdString expandedPath;
    bool isExpanded = false, wasPathSeparator = true;
    const size_t nChars = pathIn.length();
    for (size_t charI = 0; charI < nChars;)
      {
      char c = pathIn[charI];
      switch (c)
        {
        case '$': // $-variable expansion
          {
          vtkStdString variable;
          while (++charI < nChars && (isalnum(pathIn[charI]) || pathIn[charI]
              == '_'))
            {
            variable += pathIn[charI];
            }
          if (variable == "FOAM_CASE") // discard path until the variable
            {
            expandedPath = this->CasePath;
            wasPathSeparator = true;
            isExpanded = true;
            }
          else if (variable == "FOAM_CASENAME")
            {
            expandedPath += this->ExtractCaseName();
            wasPathSeparator = false;
            }
          else
            {
            const char *value = getenv(variable.c_str());
            if (value != NULL)
              {
              expandedPath += value;
              }
            const vtkStdString::size_type len = expandedPath.length();
            if (len > 0)
              {
              const char c2 = expandedPath[len - 1];
              wasPathSeparator = (c2 == '/' || c2 == '\\');
              }
            else
              {
              wasPathSeparator = false;
              }
            }
          }
          break;
        case '~': // home directory expansion
          // not using vtksys::SystemTools::ConvertToUnixSlashes() for
          // a bit better handling of "~"
          if (wasPathSeparator)
            {
            vtkStdString userName;
            while (++charI < nChars && (pathIn[charI] != '/' && pathIn[charI]
                != '\\') && pathIn[charI] != '$')
              {
              userName += pathIn[charI];
              }
            if (userName == "")
              {
              const char *homePtr = getenv("HOME");
              if (homePtr == NULL)
                {
#if defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__)
                expandedPath = "";
#else
                const struct passwd *pwentry = getpwuid(getuid());
                if (pwentry == NULL)
                  {
                  throw this->StackString() << "Home directory path not found";
                  }
                expandedPath = pwentry->pw_dir;
#endif
                }
              else
                {
                expandedPath = homePtr;
                }
              }
            else
              {
#if defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__)
              const char *homePtr = getenv("HOME");
              expandedPath
              = this->ExtractPath(homePtr ? homePtr : "") + userName;
#else
              if (userName == "OpenFOAM")
                {
                // so far only "~/.OpenFOAM" expansion is supported
                const char *homePtr = getenv("HOME");
                if (homePtr == NULL)
                  {
                  expandedPath = "";
                  }
                else
                  {
                  expandedPath = vtkStdString(homePtr) + "/.OpenFOAM";
                  }
                }
              else
                {
                const struct passwd *pwentry = getpwnam(userName.c_str());
                if (pwentry == NULL)
                  {
                  throw this->StackString() << "Home directory for user "
                  << userName.c_str() << " not found";
                  }
                expandedPath = pwentry->pw_dir;
                }
#endif
              }
            wasPathSeparator = false;
            isExpanded = true;
            break;
            }
          // fall through
        default:
          wasPathSeparator = (c == '/' || c == '\\');
          expandedPath += c;
          charI++;
        }
      }
    if (isExpanded || expandedPath.substr(0, 1) == "/" || expandedPath.substr(
        0, 1) == "\\")
      {
      return expandedPath;
      }
    else
      {
      return defaultPath + expandedPath;
      }
  }

  void IncludeFile(const vtkStdString& includedFileName,
    const vtkStdString& defaultPath, const bool includeIfPresent)
  {
    if (this->StackI >= VTK_FOAMFILE_INCLUDE_STACK_SIZE)
      {
      throw this->StackString() << "Exceeded maximum #include recursions of "
      << VTK_FOAMFILE_INCLUDE_STACK_SIZE;
      }
    // use the default bitwise copy constructor
    this->Stack[this->StackI++] = new vtkFoamFileStack(*this);
    this->Superclass::Reset();

    if(includeIfPresent)
      {
      // the behavior of #includeIfPresent is in fact "includeIfOpenable"
      try
        {
        this->Open(this->ExpandPath(includedFileName, defaultPath));
        }
      catch(vtkFoamError)
        {
        this->CloseIncludedFile();
        }
      }
    else
      {
      this->Open(this->ExpandPath(includedFileName, defaultPath));
      }
  }

  // the tokenizer
  // returns true if success, false if encountered EOF
  bool Read(vtkFoamToken& token)
  {
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == 47) // '/' == 47
      {
      this->PutBack(c);
      c = this->NextTokenHead();
      }
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    if(c != '#')
      {
      this->Superclass::WasNewline = false;
      }
#endif

    const int MAXLEN = 1024;
    char buf[MAXLEN + 1];
    int charI = 0;
    switch (c)
      {
      case '(':
      case ')':
        // high-priority punctuation token
        token = static_cast<char>(c);
        return true;
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
      case '0':
      case '-':
        // undetermined number token
        do
          {
          buf[charI++] = c;
          } while (isdigit(c = this->Getc()) && charI < MAXLEN);
        if (c != '.' && c != 'e' && c != 'E' && charI < MAXLEN && c != EOF)
          {
          // label token
          buf[charI] = '\0';
          token = static_cast<int>(strtol(buf, NULL, 10));
          this->PutBack(c);
          return true;
          }
        // fall through
      case '.':
        // scalar token
        if (c == '.' && charI < MAXLEN)
          {
          // read decimal fraction part
          buf[charI++] = c;
          while (isdigit(c = this->Getc()) && charI < MAXLEN)
            {
            buf[charI++] = c;
            }
          }
        if ((c == 'e' || c == 'E') && charI < MAXLEN)
          {
          // read exponent part
          buf[charI++] = c;
          if (((c = this->Getc()) == '+' || c == '-') && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = this->Getc();
            }
          while (isdigit(c) && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = this->Getc();
            }
          }
        if (charI == 1 && buf[0] == '-')
          {
          token = '-';
          this->PutBack(c);
          return true;
          }
        buf[charI] = '\0';
#if VTK_FOAMFILE_LOCALE_WORKAROUND
          {
          vtksys_ios::istringstream conversionStream(buf);
          double value;
          conversionStream >> value;
          token = value;
          }
#else
        token = strtod(buf, NULL);
#endif
        this->PutBack(c);
        break;
      case ';':
      case '{':
      case '}':
      case '[':
      case ']':
      case ':':
      case ',':
      case '=':
      case '+':
      case '*':
      case '/':
        // low-priority punctuation token
        token = static_cast<char>(c);
        return true;
      case '"':
        {
        // string token
        bool wasEscape = false;
        while ((c = this->Getc()) != EOF && charI < MAXLEN)
          {
          if (c == '\\' && !wasEscape)
            {
            wasEscape = true;
            continue;
            }
          else if (c == '"' && !wasEscape)
            {
            break;
            }
          else if (c == '\n')
            {
            ++this->Superclass::LineNumber;
            if (!wasEscape)
              {
              throw this->StackString()
              << "Unescaped newline in string constant";
              }
            }
          else if (wasEscape && c != '"' && c != '\n')
            {
            buf[charI++] = '\\';
            if (charI >= MAXLEN)
              {
              break;
              }
            }
          buf[charI++] = c;
          wasEscape = false;
          }
        buf[charI] = '\0';
        token = buf;
        }
        break;
      case EOF:
        // end of file
        token.SetBad();
        return false;
      case '$':
        {
        vtkFoamToken identifierToken;
        if (!this->Read(identifierToken))
          {
          throw this->StackString() << "Unexpected EOF reading identifier";
          }
        if (identifierToken.GetType() != vtkFoamToken::WORD)
          {
          throw this->StackString() << "Expected a word, found "
          << identifierToken;
          }
        token.SetIdentifier(identifierToken.ToStdString());
        return true;
        }
      case '#':
        {
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        // placing #-directives in the middle of a line looks like
        // valid for the genuine OF 1.5 parser
        if(!this->Superclass::WasNewline)
          {
          throw this->StackString()
          << "Encountered #-directive in the middle of a line";
          }
        this->Superclass::WasNewline = false;
#endif
        // read directive
        vtkFoamToken directiveToken;
        if (!this->Read(directiveToken))
          {
          throw this->StackString() << "Unexpected EOF reading directive";
          }
        if (directiveToken == "include" || directiveToken == "includeIfPresent")
          {
          vtkFoamToken fileNameToken;
          if (!this->Read(fileNameToken))
            {
            throw this->StackString() << "Unexpected EOF reading filename";
            }
          if(fileNameToken.GetType() != vtkFoamToken::STRING)
            {
            throw this->StackString() << "Expected string, found " << fileNameToken;
            }
          this->IncludeFile(fileNameToken.ToStdString(), this->GetFilePath(),
            directiveToken == "includeIfPresent");
          }
        else if (directiveToken == "inputMode")
          {
          vtkFoamToken modeToken;
          if (!this->Read(modeToken))
            {
            throw this->StackString()
            << "Unexpected EOF reading inputMode specifier";
            }
          if (modeToken == "merge" || modeToken == "default")
            {
            this->InputMode = INPUT_MODE_MERGE;
            }
          else if (modeToken == "protect")
            {
            this->InputMode = INPUT_MODE_PROTECT;
            }
          else if (modeToken == "overwrite")
            {
            this->InputMode = INPUT_MODE_OVERWRITE;
            }
          else if (modeToken == "warn")
            {
            this->InputMode = INPUT_MODE_WARN;
            }
          else if (modeToken == "error")
            {
            this->InputMode = INPUT_MODE_ERROR;
            }
          else
            {
            throw this->StackString() << "Expected one of inputMode specifiers "
            "(merge, protect, overwrite, warn, error, default), found " << modeToken;
            }
          }
        else
          {
          throw this->StackString() << "Unsupported directive "
          << directiveToken;
          }
        return this->Read(token);
        }
      default:
        // word token
        int inBrace = 0;
        do
          {
          if (c == '(')
            {
            inBrace++;
            }
          else if (c == ')' && --inBrace == -1)
            {
            break;
            }
          buf[charI++] = c;
          // valid characters that constitutes a word
          // cf. src/OpenFOAM/primitives/strings/word/wordI.H
          } while ((c = this->Getc()) != EOF && !isspace(c) && c != '"' && c
            != '/' && c != ';' && c != '{' && c != '}' && charI < MAXLEN);
        buf[charI] = '\0';
        token.SetWord(buf);
        this->PutBack(c);
      }

    if (c == EOF)
      {
      this->ThrowUnexpectedEOFException();
      }
    if (charI == MAXLEN)
      {
      throw this->StackString() << "Exceeded maximum allowed length of "
      << MAXLEN << " chars";
      }
    return true;
  }

  void Open(const vtkStdString& fileName)
  {
    // reset line number to indicate the beginning of the file when an
    // exception is thrown
    this->Superclass::LineNumber = 0;
    this->Superclass::FileName = fileName;

    if (this->Superclass::File)
      {
      throw this->StackString() << "File already opened within this object";
      }

    if ((this->Superclass::File = fopen(this->Superclass::FileName.c_str(),
        "rb")) == NULL)
      {
      throw this->StackString() << "Can't open";
      }

    unsigned char zMagic[2];
    if (fread(zMagic, 1, 2, this->Superclass::File) == 2 && zMagic[0] == 0x1f
        && zMagic[1] == 0x8b)
      {
      // gzip-compressed format
      this->Superclass::Z.avail_in = 0;
      this->Superclass::Z.next_in = Z_NULL;
      // + 32 to automatically recognize gzip format
      if (inflateInit2(&this->Superclass::Z, 15 + 32) == Z_OK)
        {
        this->Superclass::IsCompressed = true;
        this->Superclass::Inbuf = new unsigned char[VTK_FOAMFILE_INBUFSIZE];
        }
      else
        {
        fclose(this->Superclass::File);
        this->Superclass::File = NULL;
        throw this->StackString() << "Can't init zstream "
        << (this->Superclass::Z.msg ? this->Superclass::Z.msg : "");
        }
      }
    else
      {
      // uncompressed format
      this->Superclass::IsCompressed = false;
      }
    rewind(this->Superclass::File);

    this->Superclass::ZStatus = Z_OK;
    this->Superclass::Outbuf = new unsigned char[VTK_FOAMFILE_OUTBUFSIZE + 1];
    this->Superclass::BufPtr = this->Superclass::Outbuf + 1;
    this->Superclass::BufEndPtr = this->Superclass::BufPtr;
    this->Superclass::LineNumber = 1;
  }

  void Close()
  {
    while (this->CloseIncludedFile())
      ;
    this->Clear();
  }

  // gzread with buffering handling
  int Read(unsigned char *buf, const int len)
  {
    int readlen;
    const int buflen = this->Superclass::BufEndPtr - this->Superclass::BufPtr;
    if (len > buflen)
      {
      memcpy(buf, this->Superclass::BufPtr, buflen);
      readlen = this->InflateNext(buf + buflen, len - buflen);
      if (readlen >= 0)
        {
        readlen += buflen;
        }
      else
        {
        if (buflen == 0) // return EOF
          {
          readlen = -1;
          }
        else
          {
          readlen = buflen;
          }
        }
      this->Superclass::BufPtr = this->Superclass::BufEndPtr;
      }
    else
      {
      memcpy(buf, this->Superclass::BufPtr, len);
      this->Superclass::BufPtr += len;
      readlen = len;
      }
    for (int i = 0; i < readlen; i++)
      {
      if (buf[i] == '\n')
        {
        this->Superclass::LineNumber++;
        }
      }
    return readlen;
  }

  void ReadExpecting(const char expected)
  {
    // skip prepending invalid chars
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == 47) // '/' == 47
      {
      this->PutBack(c);
      c = this->NextTokenHead();
      }
    if (c != expected)
      {
      this->ThrowUnexpectedTokenException(expected, c);
      }
  }

  void ReadExpecting(const char* str)
  {
    vtkFoamToken t;
    if (!this->Read(t) || t != str)
      {
      throw this->StackString() << "Expected word \"" << str << "\", found "
      << t;
      }
  }

  int ReadIntValue();
  float ReadFloatValue();
};

int vtkFoamFile::ReadNext()
{
  if (!this->InflateNext(this->Superclass::Outbuf + 1, VTK_FOAMFILE_OUTBUFSIZE))
    {
    return this->CloseIncludedFile() ? this->Getc() : EOF;
    }
  return *this->Superclass::BufPtr++;
}

// specialized for reading an integer value.
// not using the standard strtol() for speed reason.
int vtkFoamFile::ReadIntValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
    {
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }
  if (c == 47) // '/' == 47
    {
    this->PutBack(c);
    c = this->NextTokenHead();
    }

  int nonNegative = c - 45; // '-' == 45
  if (nonNegative == 0 || c == 43) // '+' == 43
    {
    c = this->Getc();
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }

  if (!isdigit(c)) // isdigit() accepts -1 as EOF
    {
    if (c == EOF)
      {
      this->ThrowUnexpectedEOFException();
      }
    else
      {
      this->PutBack(c);
      vtkFoamToken token;
      this->Read(token);
      this->ThrowUnexpectedNonNumberTokenException(token);
      }
    }

  int num = c - 48; // '0' == 48
  while (isdigit(c = this->Getc()))
    {
    num = 10 * num + c - 48;
    }

  if (c == EOF)
    {
    this->ThrowUnexpectedEOFException();
    }
  this->PutBack(c);

  return nonNegative ? num : -num;
}

// extreamely simplified high-performing string to floating point
// conversion code based on
// ParaView3/VTK/Utilities/vtksqlite/vtk_sqlite3.c
float vtkFoamFile::ReadFloatValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
    {
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }
  if (c == 47) // '/' == 47
    {
    this->PutBack(c);
    c = this->NextTokenHead();
    }

  // determine sign
  int nonNegative = c - 45; // '-' == 45
  if (nonNegative == 0 || c == 43) // '+' == 43
    {
    c = this->Getc();
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }

  if (!isdigit(c) && c != 46) // '.' == 46, isdigit() accepts EOF
    {
    this->PutBack(c);
    vtkFoamToken token;
    this->Read(token);
    if(token == "nan" || token == "NaN")
      {
      return static_cast<float>(vtkMath::Nan());
      }
    else if(token == "inf" || token == "Inf")
      {
      return static_cast<float>(nonNegative ? vtkMath::Inf() : vtkMath::NegInf());
      }
    this->ThrowUnexpectedNonNumberTokenException(token);
    }

  // read integer part
  double num = c - 48; // '0' == 48
  while (isdigit(c = this->Getc()))
    {
    num = num * 10.0 + (c - 48);
    }

  // read decimal part
  if (c == 46) // '.'
    {
    double divisor = 1.0;

    while (isdigit(c = this->Getc()))
      {
      num = num * 10.0 + (c - 48);
      divisor *= 10.0;
      }
    num /= divisor;
    }

  // read exponent part
  if (c == 69 || c == 101) // 'E' == 69, 'e' == 101
    {
    int esign = 1;
    int eval = 0;
    double scale = 1.0;

    c = this->Getc();
    if (c == 45) // '-'
      {
      esign = -1;
      c = this->Getc();
      }
    else if (c == 43) // '+'
      {
      c = this->Getc();
      }

    while (isdigit(c))
      {
      eval = eval * 10 + (c - 48);
      c = this->Getc();
      }

    // fast exponent multiplication!
    while (eval >= 64)
      {
      scale *= 1.0e+64;
      eval -= 64;
      }
    while (eval >= 16)
      {
      scale *= 1.0e+16;
      eval -= 16;
      }
    while (eval >= 4)
      {
      scale *= 1.0e+4;
      eval -= 4;
      }
    while (eval >= 1)
      {
      scale *= 1.0e+1;
      eval -= 1;
      }

    if (esign < 0)
      {
      num /= scale;
      }
    else
      {
      num *= scale;
      }
    }

  if (c == EOF)
    {
    this->ThrowUnexpectedEOFException();
    }
  this->PutBack(c);

  return static_cast<float>(nonNegative ? num : -num);
}

// hacks to keep exception throwing code out-of-line to make
// putBack() and readExpecting() inline expandable
void vtkFoamFile::ThrowUnexpectedEOFException()
{
  throw this->StackString() << "Unexpected EOF";
}

void vtkFoamFile::ThrowUnexpectedNonNumberTokenException(const vtkFoamToken &token)
{
  throw this->StackString() << "Expected a number, found a non-number token "
  << token;
}

void vtkFoamFile::ThrowUnexpectedTokenException(const char expected, const int c)
{
  vtkFoamError sstr;
  sstr << this->StackString() << "Expected punctuation token '" << expected
      << "', found ";
  if (c == EOF)
    {
    sstr << "EOF";
    }
  else
    {
    sstr << static_cast<char>(c);
    }
  throw sstr;
}

void vtkFoamFile::ThrowDuplicatedPutBackException()
{
  throw this->StackString() << "Attempted duplicated putBack()";
}

bool vtkFoamFile::InflateNext(unsigned char *buf,
    int requestSize)
{
  size_t size;
  if (this->Superclass::IsCompressed)
    {
    if (this->Superclass::ZStatus != Z_OK)
      {
      return false;
      }
    this->Superclass::Z.next_out = buf;
    this->Superclass::Z.avail_out = requestSize;

    do
      {
      if (this->Superclass::Z.avail_in == 0)
        {
        this->Superclass::Z.next_in = this->Superclass::Inbuf;
        this->Superclass::Z.avail_in = static_cast<uInt>(fread(this->Superclass::Inbuf, 1,
            VTK_FOAMFILE_INBUFSIZE, this->Superclass::File));
        if (ferror(this->Superclass::File))
          {
          throw this->StackString() << "Fread failed";
          }
        }
      this->Superclass::ZStatus = inflate(&this->Superclass::Z, Z_NO_FLUSH);
      if (this->Superclass::ZStatus == Z_STREAM_END
#if VTK_FOAMFILE_OMIT_CRCCHECK
      // the dummy CRC function causes data error when finalizing
      // so we have to proceed even when a data error is detected
      || this->Superclass::ZStatus == Z_DATA_ERROR
#endif
      )
        {
        break;
        }
      if (this->Superclass::ZStatus != Z_OK)
        {
        throw this->StackString() << "Inflation failed: "
        << (this->Superclass::Z.msg ? this->Superclass::Z.msg : "");
        }
      } while (this->Superclass::Z.avail_out > 0);

    size = requestSize - this->Superclass::Z.avail_out;
    }
  else
    {
    // not compressed
    size = fread(buf, 1, requestSize, this->Superclass::File);
    }

  if (size <= 0)
    {
    // retain the current location bufPtr_ to the end of the buffer so that
    // getc() returns EOF again when called next time
    return false;
    }
  // size > 0
  // reserve the first byte for getback char
  this->Superclass::BufPtr = this->Superclass::Outbuf + 1;
  this->Superclass::BufEndPtr = this->Superclass::BufPtr + size;
  return true;
}

// get next semantically valid character
int vtkFoamFile::NextTokenHead()
{
  for (;;)
    {
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == '/')
      {
      if ((c = this->Getc()) == '/')
        {
        while ((c = this->Getc()) != EOF && c != '\n')
          ;
        if (c == EOF)
          {
          return c;
          }
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      else if (c == '*')
        {
        for (;;)
          {
          while ((c = this->Getc()) != EOF && c != '*')
            {
            if (c == '\n')
              {
              ++this->Superclass::LineNumber;
              }
            }
          if (c == EOF)
            {
            return c;
            }
          else if ((c = this->Getc()) == '/')
            {
            break;
            }
          this->PutBack(c);
          }
        }
      else
        {
        this->PutBack(c); // may be an EOF
        return '/';
        }
      }
    else // may be an EOF
      {
      return c;
      }
    }
#if defined(__hpux)
  return EOF; // this line should not be executed; workaround for HP-UXia64-aCC
#endif
}

//-----------------------------------------------------------------------------
// class vtkFoamIOobject
// holds file handle, file format, name of the object the file holds and
// type of the object.
struct vtkFoamIOobject : public vtkFoamFile
{
private:
  typedef vtkFoamFile Superclass;

public:
  enum fileFormat
    {UNDEFINED, ASCII, BINARY};

private:
  fileFormat Format;
  vtkStdString ObjectName;
  vtkStdString HeaderClassName;
  bool Is13Positions;
  bool IsSinglePrecisionBinary;
  vtkFoamError E;

  vtkFoamIOobject();
  void ReadHeader(); // defined later
public:
  vtkFoamIOobject(const vtkStdString& casePath, const bool isSinglePrecisionBinary) :
    vtkFoamFile(casePath), Format(UNDEFINED), Is13Positions(false),
    IsSinglePrecisionBinary(isSinglePrecisionBinary), E()
  {
  }
  ~vtkFoamIOobject()
  {
    this->Close();
  }

  bool Open(const vtkStdString& file)
  {
    try
      {
      this->Superclass::Open(file);
      }
    catch(vtkFoamError& e)
      {
      this->E = e;
      return false;
      }

    try
      {
      this->ReadHeader();
      }
    catch(vtkFoamError& e)
      {
      this->Superclass::Close();
      this->E = e;
      return false;
      }
    return true;
  }

  void Close()
  {
    this->Superclass::Close();
    this->Format = UNDEFINED;
    this->ObjectName.erase();
    this->HeaderClassName.erase();
    this->E.erase();
  }
  fileFormat GetFormat() const
  {
    return this->Format;
  }
  const vtkStdString& GetClassName() const
  {
    return this->HeaderClassName;
  }
  const vtkStdString& GetObjectName() const
  {
    return this->ObjectName;
  }
  const vtkFoamError& GetError() const
  {
    return this->E;
  }
  void SetError(const vtkFoamError& e)
  {
    this->E = e;
  }
  void SetIs13Positions(const bool is13Positions)
  {
    this->Is13Positions = is13Positions;
  }
  bool GetIs13Positions() const
  {
    return this->Is13Positions;
  }
  bool GetIsSinglePrecisionBinary() const
  {
    return this->IsSinglePrecisionBinary;
  }
};

//-----------------------------------------------------------------------------
// workarounding class for older compilers (gcc-3.3.x and possibly older)
template <typename T> struct vtkFoamReadValue
{
public:
  static T ReadValue(vtkFoamIOobject &io);
};

VTK_TEMPLATE_SPECIALIZE inline int vtkFoamReadValue<int>::ReadValue(vtkFoamIOobject& io)
{
  return io.ReadIntValue();
}

VTK_TEMPLATE_SPECIALIZE inline float vtkFoamReadValue<float>::ReadValue(vtkFoamIOobject& io)
{
  return io.ReadFloatValue();
}

//-----------------------------------------------------------------------------
// class vtkFoamEntryValue
// a class that represents a value of a dictionary entry that corresponds to
// its keyword. note that an entry can have more than one value.
struct vtkFoamEntryValue : public vtkFoamToken
{
public:
  enum uniformTypes
    {UNDEFINED, UNIFORM, NONUNIFORM};
private:
  typedef vtkFoamToken Superclass;

  uniformTypes IsUniform;
  bool Managed;
  const vtkFoamEntry *UpperEntryPtr;

  vtkFoamEntryValue();
  vtkObjectBase *ToVTKObject()
  {
    return this->Superclass::VtkObjectPtr;
  }
  void Clear();
  void ReadList(vtkFoamIOobject& io);

public:
  // reads primitive int/float lists
  template <typename listT, typename primitiveT> class listTraits
  {
    listT *Ptr;

  public:
    listTraits() :
      Ptr(listT::New())
    {
    }
    listT *GetPtr()
    {
      return this->Ptr;
    }
    void ReadUniformValues(vtkFoamIOobject& io, const int size)
    {
      primitiveT value = vtkFoamReadValue<primitiveT>::ReadValue(io);
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetValue(i, value);
        }
    }
    void ReadAsciiList(vtkFoamIOobject& io, const int size)
    {
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetValue(i, vtkFoamReadValue<primitiveT>::ReadValue(io));
        }
    }
    void ReadBinaryList(vtkFoamIOobject& io, const int size)
    {
      io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)), size
          * sizeof(primitiveT));
    }
    void ReadValue(vtkFoamIOobject&, vtkFoamToken& currToken)
    {
      if (!currToken.Is<primitiveT>())
        {
        throw vtkFoamError() << "Expected an integer or a (, found "
        << currToken;
        }
      this->Ptr->InsertNextValue(currToken.To<primitiveT>());
    }
  };

  // reads rank 1 lists of types vector, sphericalTensor, symmTensor
  // and tensor. if isPositions is true it reads Cloud type of data as
  // particle positions. cf. (the positions format)
  // src/lagrangian/basic/particle/particleIO.C
  template <typename listT, typename primitiveT, int nComponents,
      bool isPositions> class vectorListTraits
  {
    listT *Ptr;

  public:
    vectorListTraits() :
      Ptr(listT::New())
    {
      this->Ptr->SetNumberOfComponents(nComponents);
    }
    listT *GetPtr()
    {
      return this->Ptr;
    }
    void ReadUniformValues(vtkFoamIOobject& io, const int size)
    {
      io.ReadExpecting('(');
      primitiveT vectorValue[nComponents];
      for (int j = 0; j < nComponents; j++)
        {
        vectorValue[j] = vtkFoamReadValue<primitiveT>::ReadValue(io);
        }
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetTuple(i, vectorValue);
        }
      io.ReadExpecting(')');
      if (isPositions)
        {
        // skip label celli
        vtkFoamReadValue<int>::ReadValue(io);
        }
    }
    void ReadAsciiList(vtkFoamIOobject& io, const int size)
    {
      for (int i = 0; i < size; i++)
        {
        io.ReadExpecting('(');
        primitiveT *vectorTupleI = this->Ptr->GetPointer(nComponents * i);
        for (int j = 0; j < nComponents; j++)
          {
          vectorTupleI[j] = vtkFoamReadValue<primitiveT>::ReadValue(io);
          }
        io.ReadExpecting(')');
        if (isPositions)
          {
          // skip label celli
          vtkFoamReadValue<int>::ReadValue(io);
          }
        }
    }
    void ReadBinaryList(vtkFoamIOobject& io, const int size)
    {
      if (isPositions) // lagrangian/positions (class Cloud)
        {
        // allocate space along with the larger 1.4 format since the
        // size must be determined at compile-time. The buffer is
        // allocated on the stack in order to avoid leak when an
        // exception is thrown.
        unsigned char buffer[sizeof(double) * (nComponents + 1) + 2
            * sizeof(int)];
        const int realSize = io.GetIsSinglePrecisionBinary()
            ? sizeof(float) : sizeof(double);
        const int nBytes = (io.GetIs13Positions()
            // skip label celli
            ? realSize * nComponents + sizeof(int)
            // skip label celli, label facei and scalar stepFraction
            : realSize * (nComponents + 1) + 2 * sizeof(int));
        for (int i = 0; i < size; i++)
          {
          io.ReadExpecting('(');
          io.Read(buffer, nBytes);
          this->Ptr->SetTuple(i, reinterpret_cast<double *>(buffer));
          io.ReadExpecting(')');
          }
        }
      else
        {
        if (io.GetIsSinglePrecisionBinary())
          {
          io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)),
              sizeof(float) * nComponents * size);
          }
        else
          {
#if 0
          // reduce calls to io.Read() by using the array space as buffer
          float *destination = this->Ptr->GetPointer(0);
          double *source = reinterpret_cast<double *>(destination);

          int remainingSize = size;
          for (int halfSize = remainingSize / 2; halfSize > 0;
               halfSize = remainingSize / 2)
            {
            const int convertSize = nComponents * halfSize;
            io.Read(reinterpret_cast<unsigned char *>(source),
            sizeof(double) * convertSize);
            for (int i = 0; i < convertSize; i++)
              {
              *destination++ = static_cast<float>(*source++);
              }
            source = reinterpret_cast<double *>(destination);
            remainingSize -= halfSize;
            }
          if (remainingSize > 0)
            {
            // read the last one element that doesn't fit into the
            // remaining space
            double buffer[nComponents];
            io.Read(reinterpret_cast<unsigned char *>(buffer),
                sizeof(double) * nComponents);
            for (int i = 0; i < nComponents; i++)
              {
              *destination++ = static_cast<float>(buffer[i]);
              }
            }
#else
          const int bufferUnit = 32, nDivs = size / bufferUnit,
            unitSize = nComponents * bufferUnit;
          double buffer[nComponents * bufferUnit];

          for (int i = 0; i < nDivs; i++)
            {
            float *destination = this->Ptr->GetPointer(i * unitSize);
            io.Read(reinterpret_cast<unsigned char *>(buffer),
                sizeof(double) * unitSize);
            for (int j = 0; j < unitSize; j++)
              {
              destination[j] = static_cast<float>(buffer[j]);
              }
            }
          const int remainingSize = nComponents * (size % bufferUnit);
          float *destination = this->Ptr->GetPointer(nDivs * unitSize);
          io.Read(reinterpret_cast<unsigned char *>(buffer),
              sizeof(double) * remainingSize);
          for (int j = 0; j < remainingSize; j++)
            {
            destination[j] = static_cast<float>(buffer[j]);
            }
#endif
          }
        }
    }
    void ReadValue(vtkFoamIOobject& io, vtkFoamToken& currToken)
    {
      if (currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      primitiveT v[nComponents];
      for (int j = 0; j < nComponents; j++)
        {
        v[j] = vtkFoamReadValue<primitiveT>::ReadValue(io);
        }
      this->Ptr->InsertNextTuple(v);
      io.ReadExpecting(')');
    }
  };

  vtkFoamEntryValue(const vtkFoamEntry *upperEntryPtr) :
    vtkFoamToken(), IsUniform(UNDEFINED), Managed(true),
        UpperEntryPtr(upperEntryPtr)
  {
  }
  vtkFoamEntryValue(vtkFoamEntryValue&, const vtkFoamEntry *);
  ~vtkFoamEntryValue()
  {
    this->Clear();
  }

  void SetEmptyList()
  {
    this->Clear();
    this->IsUniform = NONUNIFORM;
    this->Superclass::Type = EMPTYLIST;
  }
  uniformTypes GetIsUniform() const
  {
    return this->IsUniform;
  }
  void SetIsUniform(const uniformTypes isUniform)
  {
    this->IsUniform = isUniform;
  }
  void Read(vtkFoamIOobject& io);
  void ReadDictionary(vtkFoamIOobject& io, const vtkFoamToken& firstKeyword);
  const vtkIntArray& LabelList() const
  {
    return *this->Superclass::LabelListPtr;
  }
  vtkIntArray& LabelList()
  {
    return *this->Superclass::LabelListPtr;
  }
  const vtkFoamIntVectorVector& LabelListList() const
  {
    return *this->Superclass::LabelListListPtr;
  }
  const vtkFloatArray& ScalarList() const
  {
    return *this->Superclass::ScalarListPtr;
  }
  vtkFloatArray& ScalarList()
  {
    return *this->Superclass::ScalarListPtr;
  }
  const vtkFloatArray& VectorList() const
  {
    return *this->Superclass::VectorListPtr;
  }
  const vtkFoamDict& Dictionary() const
  {
    return *this->Superclass::DictPtr;
  }
  vtkFoamDict& Dictionary()
  {
    return *this->Superclass::DictPtr;
  }

  void *Ptr()
  {
    this->Managed = false; // returned pointer will not be deleted by the d'tor
    // all list pointers are in a single union
    return (void *)this->Superclass::LabelListPtr;
  }

  vtkStdString ToStdString() const
  {
    return this->Superclass::IsWordOrString()
        || this->Superclass::GetType() == vtkFoamToken::IDENTIFIER
        ? this->Superclass::ToStdString() : vtkStdString();
  }
  vtkStdString ToWord() const
  {
    return this->Superclass::Type == WORD ? this->Superclass::ToStdString()
        : vtkStdString();
  }
  vtkStdString ToString() const
  {
    return this->Superclass::Type == STRING ? this->Superclass::ToStdString()
        : vtkStdString();
  }
  float ToFloat() const
  {
    return this->Superclass::Type == SCALAR || this->Superclass::Type == LABEL ? this->Superclass::To<float>()
        : 0.0F;
  }
  double ToDouble() const
  {
    return this->Superclass::Type == SCALAR || this->Superclass::Type == LABEL ? this->Superclass::To<double>()
        : 0.0;
  }
  int ToInt() const
  {
    return this->Superclass::Type == LABEL ? this->Superclass::To<int>() : 0;
  }

  // the following two are for an exceptional expression of
  // `LABEL{LABELorSCALAR}' without type prefix (e. g. `2{-0}' in
  // mixedRhoE B.C. in rhopSonicFoam/shockTube)
  void MakeLabelList(const int labelValue, const int size)
  {
    this->Superclass::LabelListPtr = vtkIntArray::New();
    this->Superclass::Type = LABELLIST;
    this->Superclass::LabelListPtr->SetNumberOfValues(size);
    for (int i = 0; i < size; i++)
      {
      this->Superclass::LabelListPtr->SetValue(i, labelValue);
      }
  }
  void MakeScalarList(const float scalarValue, const int size)
  {
    this->Superclass::ScalarListPtr = vtkFloatArray::New();
    this->Superclass::Type = SCALARLIST;
    this->Superclass::ScalarListPtr->SetNumberOfValues(size);
    for (int i = 0; i < size; i++)
      {
      this->Superclass::ScalarListPtr->SetValue(i, scalarValue);
      }
  }

  // reads dimensionSet
  void ReadDimensionSet(vtkFoamIOobject& io)
  {
    const int nDims = 7;
    this->Superclass::LabelListPtr = vtkIntArray::New();
    this->Superclass::Type = LABELLIST;
    this->Superclass::LabelListPtr->SetNumberOfValues(nDims);
    for (int dimI = 0; dimI < nDims; dimI++)
      {
      this->Superclass::LabelListPtr->SetValue(dimI, vtkFoamReadValue<int>::ReadValue(io));
      }
    io.ReadExpecting(']');
  }

  template <vtkFoamToken::tokenType listType, typename traitsT> void ReadNonuniformList(
      vtkFoamIOobject& io);

  // reads a list of labelLists. requires size prefix of the listList
  // to be present. size of each sublist must also be present in the
  // stream if the format is binary.
  void ReadLabelListList(vtkFoamIOobject& io)
  {
    vtkFoamToken currToken;
    if (!io.Read(currToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if (currToken.GetType() == vtkFoamToken::LABEL)
      {
      const int sizeI = currToken.To<int>();
      if (sizeI < 0)
        {
        throw vtkFoamError() << "List size must not be negative: size = "
        << sizeI;
        }
      // gives initial guess for list size. Assumes 4.25 vertices per
      // face in average, taking into account typical snappyHexMesh
      // polyhedral meshes. The +0.25 vertices per face over a quad
      // should be enough for most split-hex polyhedral meshes - even
      // the iglooWithFridges case, which appears to be a heavy
      // polyhedral mesh, only has 4.124 vertices per face. The number
      // will not be enough for polyDualMesh-type meshes though.
      this->Superclass::LabelListListPtr = new vtkFoamIntVectorVector(sizeI, 4 * sizeI + (sizeI + 3) / 4);
      this->Superclass::Type = LABELLISTLIST;
      io.ReadExpecting('(');
      int bodyI = 0;
      for (int i = 0; i < sizeI; i++)
        {
        if (!io.Read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }
        if (currToken.GetType() == vtkFoamToken::LABEL)
          {
          const int sizeJ = currToken.To<int>();
          if (sizeJ < 0)
            {
            throw vtkFoamError() << "List size must not be negative: size = "
            << sizeJ;
            }
          int *listI = this->Superclass::LabelListListPtr->WritePointer(i,
              bodyI, sizeJ);

          if (io.GetFormat() == vtkFoamIOobject::ASCII)
            {
            io.ReadExpecting('(');
            for (int j = 0; j < sizeJ; j++)
              {
              listI[j] = vtkFoamReadValue<int>::ReadValue(io);
              }
            io.ReadExpecting(')');
            }
          else
            {
            if (sizeJ > 0) // avoid invalid reference to labelListI.at(0)
              {
              io.ReadExpecting('(');
              io.Read(reinterpret_cast<unsigned char*>(listI), sizeJ
                  * sizeof(int));
              io.ReadExpecting(')');
              }
            }
          bodyI += sizeJ;
          }
        else if (currToken == '(')
          {
          this->Superclass::LabelListListPtr->SetIndex(i, bodyI);
          while (io.Read(currToken) && currToken != ')')
            {
            if (currToken.GetType() != vtkFoamToken::LABEL)
              {
              throw vtkFoamError() << "Expected an integer, found "
              << currToken;
              }
            this->Superclass::LabelListListPtr
                ->InsertValue(bodyI++, currToken.To<int>());
            }
          }
        else
          {
          throw vtkFoamError() << "Expected integer or '(', found "
          << currToken;
          }
        }
      // set the next index of the last element required to calculate
      // the last subarray size
      this->Superclass::LabelListListPtr->SetIndex(sizeI, bodyI);
      // shrink to the actually used size
      this->Superclass::LabelListListPtr->ResizeBody(bodyI);
      io.ReadExpecting(')');
      }
    else
      {
      throw vtkFoamError() << "Expected integer, found " << currToken;
      }
  }

  // reads compact list of labels.
  void ReadCompactIOLabelList(vtkFoamIOobject& io)
  {
    if (io.GetFormat() != vtkFoamIOobject::BINARY)
      {
      this->ReadLabelListList(io);
      return;
      }

    this->Superclass::LabelListListPtr = new vtkFoamIntVectorVector;
    this->Superclass::Type = LABELLISTLIST;
    for(int arrayI = 0; arrayI < 2; arrayI++)
      {
      vtkFoamToken currToken;
      if (!io.Read(currToken))
        {
        throw vtkFoamError() << "Unexpected EOF";
        }
      if (currToken.GetType() == vtkFoamToken::LABEL)
        {
        const int sizeI = currToken.To<int>();
        if (sizeI < 0)
          {
          throw vtkFoamError() << "List size must not be negative: size = "
              << sizeI;
          }
        if (sizeI > 0) // avoid invalid reference
          {
          vtkIntArray *array = (arrayI == 0
              ? this->Superclass::LabelListListPtr->GetIndices()
              : this->Superclass::LabelListListPtr->GetBody());
          array->SetNumberOfValues(sizeI);
          io.ReadExpecting('(');
          io.Read(reinterpret_cast<unsigned char*>(array->GetPointer(0)),
              sizeI * sizeof(int));
          io.ReadExpecting(')');
          }
        }
      else
        {
        throw vtkFoamError() << "Expected integer, found " << currToken;
        }
      }
  }

  bool ReadField(vtkFoamIOobject& io)
  {
    try
      {
      // lagrangian labels (cf. gnemdFoam/nanoNozzle)
      if(io.GetClassName() == "labelField")
        {
        this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
        }
      // lagrangian scalars

      else if(io.GetClassName() == "scalarField")
        {
        this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(
            io);
        }
      else if(io.GetClassName() == "sphericalTensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 1, false> >(io);
        }
      // polyMesh/points, lagrangian vectors

      else if(io.GetClassName() == "vectorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 3, false> >(io);
        }
      else if(io.GetClassName() == "symmTensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 6, false> >(io);
        }
      else if(io.GetClassName() == "tensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 9, false> >(io);
        }
      else
        {
        throw vtkFoamError() << "Non-supported field type "
        << io.GetClassName();
        }
      }
    catch(vtkFoamError& e)
      {
      io.SetError(e);
      return false;
      }
    return true;
  }
};

// specialization for reading double precision binary into vtkFloatArray.
// Must precede ReadNonuniformList() below (HP-UXia64-aCC).
VTK_TEMPLATE_SPECIALIZE
void vtkFoamEntryValue::listTraits<vtkFloatArray, float>::ReadBinaryList(
    vtkFoamIOobject& io, const int size)
{
  if(io.GetIsSinglePrecisionBinary())
    {
    io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)),
        sizeof(float) * size);
    }
  else
    {
#if 0
    // reduce calls to io.Read() by using the array space as buffer
    float *destination = this->Ptr->GetPointer(0);
    double *source = reinterpret_cast<double *>(destination);

    int remainingSize = size;
    for (int halfSize = size / 2; halfSize > 0; halfSize = remainingSize / 2)
      {
      io.Read(reinterpret_cast<unsigned char *>(source),
          sizeof(double) * halfSize);
      for (int i = 0; i < halfSize; i++)
        {
        *destination++ = static_cast<float>(*source++);
        }
      source = reinterpret_cast<double *>(destination);
      remainingSize -= halfSize;
      }
    if (remainingSize > 0)
      {
      // read the last one element that doesn't fit into the remaining space
      double buffer;
      io.Read(reinterpret_cast<unsigned char *>(&buffer), sizeof(double));
      *destination = static_cast<float>(buffer);
      }
#else
    const int bufferUnit = 32, nDivs = size / bufferUnit;
    double buffer[bufferUnit];

    for (int i = 0; i < nDivs; i++)
      {
      float *destination = this->Ptr->GetPointer(i * bufferUnit);
      io.Read(reinterpret_cast<unsigned char *>(buffer),
          sizeof(double) * bufferUnit);
      for (int j = 0; j < bufferUnit; j++)
        {
        destination[j] = static_cast<float>(buffer[j]);
        }
      }
    const int remainingSize = size % bufferUnit;
    float *destination = this->Ptr->GetPointer(nDivs * bufferUnit);
    io.Read(reinterpret_cast<unsigned char *>(buffer),
        sizeof(double) * remainingSize);
    for (int j = 0; j < remainingSize; j++)
      {
      destination[j] = static_cast<float>(buffer[j]);
      }
#endif
    }
}

// generic reader for nonuniform lists. requires size prefix of the
// list to be present in the stream if the format is binary.
template <vtkFoamToken::tokenType listType, typename traitsT>
void vtkFoamEntryValue::ReadNonuniformList(vtkFoamIOobject& io)
{
  vtkFoamToken currToken;
  if (!io.Read(currToken))
    {
    throw vtkFoamError() << "Unexpected EOF";
    }
  traitsT list;
  this->Superclass::Type = listType;
  this->Superclass::VtkObjectPtr = list.GetPtr();
  if (currToken.Is<int>())
    {
    const int size = currToken.To<int>();
    if (size < 0)
      {
      throw vtkFoamError() << "List size must not be negative: size = " << size;
      }
    list.GetPtr()->SetNumberOfTuples(size);
    if (io.GetFormat() == vtkFoamIOobject::ASCII)
      {
      if (!io.Read(currToken))
        {
        throw vtkFoamError() << "Unexpected EOF";
        }
      // some objects have lists with only one element enclosed by {}
      // e. g. simpleFoam/pitzDaily3Blocks/constant/polyMesh/faceZones
      if (currToken == '{')
        {
        list.ReadUniformValues(io, size);
        io.ReadExpecting('}');
        return;
        }
      else if (currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      list.ReadAsciiList(io, size);
      io.ReadExpecting(')');
      }
    else
      {
      if (size > 0)
        {
        // read parentheses only when size > 0
        io.ReadExpecting('(');
        list.ReadBinaryList(io, size);
        io.ReadExpecting(')');
        }
      }
    }
  else if (currToken == '(')
    {
    while (io.Read(currToken) && currToken != ')')
      {
      list.ReadValue(io, currToken);
      }
    list.GetPtr()->Squeeze();
    }
  else
    {
    throw vtkFoamError() << "Expected integer or '(', found " << currToken;
    }
}

//-----------------------------------------------------------------------------
// class vtkFoamKeyword
// a class that handles regex
#ifdef VTK_FOAMFILE_HAVE_REGEX
// use the system POSIX regex if it was found
struct vtkFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  regex_t *Preg;
  void operator=(const vtkFoamKeyword &); // not implemented

  void Clear()
  {
    if (this->Preg != NULL)
      {
      regfree(this->Preg);
      delete this->Preg;
      this->Preg = NULL;
      }
  }

  int Compile(const vtkStdString &keyword)
  {
    this->Preg = new regex_t;
    return regcomp(this->Preg, keyword.c_str(), REG_EXTENDED);
  }

public:
  vtkFoamKeyword()
    : vtkStdString(), Preg(NULL)
  {
  }
  vtkFoamKeyword(const vtkFoamKeyword &keyword)
    : vtkStdString(keyword), Preg(NULL)
  {
    if (keyword.Preg != NULL)
      {
      // does not handle error assuming the compilation always succeed
      // since exception should not be thrown from within the constructor
      this->Compile(keyword);
      }
  }
  ~vtkFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkFoamToken::STRING && isRegEx)
      {
      const int ret = this->Compile(keywordToken.ToStdString());
      if (ret != 0)
        {
        const int msgSize = regerror(ret, this->Preg, NULL, 0);
        char *msgBuf = new char[msgSize + 1];
        regerror(ret, this->Preg, msgBuf, msgSize + 1);
        vtkStdString errString = msgBuf;
        delete [] msgBuf;
        this->Clear();
        throw vtkFoamError() << "regular expression error: " << errString;
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Preg == NULL)
      {
      return false;
      }

    regmatch_t pmatch;
    const int ret = regexec(this->Preg, string.c_str(), 1, &pmatch, 0);
    return ret == 0 && pmatch.rm_so == 0
        && static_cast<size_t>(pmatch.rm_eo) == string.length();
  }
};
#elif defined(_MSC_VER) && (_MSC_FULL_VER > 150030729 \
    || _MSC_FULL_VER == 150030729 && _MSC_BUILD >= 1)
// if VS2008SP1 or above, use std::tr1::regex
#include <regex>
struct vtkFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  vtkstd::tr1::regex *Regex;
  void operator=(const vtkFoamKeyword &); // not implemented

  void Clear()
  {
    delete this->Regex;
    this->Regex = NULL;
  }

public:
  vtkFoamKeyword()
    : vtkStdString(), Regex(NULL)
  {
  }
  vtkFoamKeyword(const vtkFoamKeyword &keyword)
    : vtkStdString(keyword), Regex(NULL)
  {
    if (keyword.Regex != NULL)
      {
      this->Regex = new vtkstd::tr1::regex(*keyword.Regex);
      }
  }
  ~vtkFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkFoamToken::STRING && isRegEx)
      {
      try
        {
        this->Regex = new vtkstd::tr1::regex(
            keywordToken.ToStdString(), vtkstd::tr1::regex::extended);
        }
      catch (const vtkstd::tr1::regex_error)
        {
        this->Clear();
        throw vtkFoamError() << "regular expression error";
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Regex == NULL)
      {
      return false;
      }

    return regex_match(string, *Regex);
  }
};
#else
// if none is found, use vtksys::RegularExpression as the fallback option
struct vtkFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  vtksys::RegularExpression *Regex;
  void operator=(const vtkFoamKeyword &); // not implemented

  void Clear()
  {
    delete this->Regex;
    this->Regex = NULL;
  }

public:
  vtkFoamKeyword()
    : vtkStdString(), Regex(NULL)
  {
  }
  vtkFoamKeyword(const vtkFoamKeyword &keyword)
    : vtkStdString(keyword), Regex(NULL)
  {
    if (keyword.Regex != NULL)
      {
      this->Regex = new vtksys::RegularExpression(*keyword.Regex);
      }
  }
  ~vtkFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkFoamToken::STRING && isRegEx)
      {
      this->Regex = new vtksys::RegularExpression(
          keywordToken.ToStdString().c_str());
      if (!this->Regex->is_valid())
        {
        this->Clear();
        throw vtkFoamError() << "regular expression error";
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Regex == NULL)
      {
      return false;
      }

    return this->Regex->find(string) && this->Regex->start() == 0
        && this->Regex->end() == string.length();
  }
};
#endif

//-----------------------------------------------------------------------------
// class vtkFoamEntry
// a class that represents an entry of a dictionary. note that an
// entry can have more than one value.
struct vtkFoamEntry : public vtkstd::vector<vtkFoamEntryValue*>
{
private:
  typedef vtkstd::vector<vtkFoamEntryValue*> Superclass;
  vtkFoamKeyword Keyword;
  vtkFoamDict *UpperDictPtr;

  vtkFoamEntry();

public:
  vtkFoamEntry(vtkFoamDict *upperDictPtr) :
    Keyword(), UpperDictPtr(upperDictPtr)
  {
  }
  vtkFoamEntry(const vtkFoamEntry& entry, vtkFoamDict *upperDictPtr) :
    Superclass(entry.size()), Keyword(entry.GetKeyword()),
        UpperDictPtr(upperDictPtr)
  {
    for (size_t valueI = 0; valueI < entry.size(); valueI++)
      {
      this->Superclass::operator[](valueI) = new vtkFoamEntryValue(*entry[valueI], this);
      }
  }

  ~vtkFoamEntry()
  {
    this->Clear();
  }

  void Clear()
  {
    for (size_t i = 0; i < this->Superclass::size(); i++)
      {
      delete this->Superclass::operator[](i);
      }
    this->Superclass::clear();
  }
  const vtkFoamKeyword& GetKeyword() const
  {
    return this->Keyword;
  }
  void SetKeyword(const vtkFoamToken& keyword, const bool isRegExKeyword)
  {
    this->Keyword.SetKeyword(keyword, isRegExKeyword);
  }
  const vtkFoamEntryValue& FirstValue() const
  {
    return *this->Superclass::operator[](0);
  }
  vtkFoamEntryValue& FirstValue()
  {
    return *this->Superclass::operator[](0);
  }
  const vtkIntArray& LabelList() const
  {
    return this->FirstValue().LabelList();
  }
  vtkIntArray& LabelList()
  {
    return this->FirstValue().LabelList();
  }
  const vtkFoamIntVectorVector& LabelListList() const
  {
    return this->FirstValue().LabelListList();
  }
  const vtkFloatArray& ScalarList() const
  {
    return this->FirstValue().ScalarList();
  }
  vtkFloatArray& ScalarList()
  {
    return this->FirstValue().ScalarList();
  }
  const vtkFloatArray& VectorList() const
  {
    return this->FirstValue().VectorList();
  }
  const vtkFoamDict& Dictionary() const
  {
    return this->FirstValue().Dictionary();
  }
  vtkFoamDict& Dictionary()
  {
    return this->FirstValue().Dictionary();
  }
  void *Ptr()
  {
    return this->FirstValue().Ptr();
  }
  const vtkFoamDict *GetUpperDictPtr() const
  {
    return this->UpperDictPtr;
  }

  vtkStdString ToStdString() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToStdString() : vtkStdString();
  }
  vtkStdString ToWord() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToWord() : vtkStdString();
  }
  vtkStdString ToString() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToString() : vtkStdString();
  }
  float ToFloat() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToFloat() : 0.0F;
  }
  double ToDouble() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToDouble() : 0.0;
  }
  int ToInt() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToInt() : 0;
  }

  void ReadDictionary(vtkFoamIOobject& io)
  {
    this->Superclass::push_back(new vtkFoamEntryValue(this));
    this->Superclass::back()->ReadDictionary(io, vtkFoamToken());
  }

  // read values of an entry
  void Read(vtkFoamIOobject& io);
};

//-----------------------------------------------------------------------------
// class vtkFoamDict
// a class that holds a FoamFile data structure
struct vtkFoamDict : public vtkstd::vector<vtkFoamEntry*>
{
private:
  typedef vtkstd::vector<vtkFoamEntry*> Superclass;

  vtkFoamToken Token;
  const vtkFoamDict *UpperDictPtr;

  vtkFoamDict(const vtkFoamDict &);

public:
  vtkFoamDict(const vtkFoamDict *upperDictPtr = NULL) :
    Superclass(), Token(), UpperDictPtr(upperDictPtr)
  {
  }
  vtkFoamDict(const vtkFoamDict& dict, const vtkFoamDict *upperDictPtr) :
    Superclass(dict.size()), Token(), UpperDictPtr(upperDictPtr)
  {
    if (dict.GetType() == vtkFoamToken::DICTIONARY)
      {
      for (size_t entryI = 0; entryI < dict.size(); entryI++)
        {
        this->operator[](entryI) = new vtkFoamEntry(*dict[entryI], this);
        }
      }
  }

  ~vtkFoamDict()
  {
    if (this->Token.GetType() == vtkFoamToken::UNDEFINED)
      {
      for (size_t i = 0; i < this->Superclass::size(); i++)
        {
        delete this->operator[](i);
        }
      }
  }

  vtkFoamToken::tokenType GetType() const
  {
    return this->Token.GetType() == vtkFoamToken::UNDEFINED ? vtkFoamToken::DICTIONARY
        : this->Token.GetType();
  }
  const vtkFoamToken &GetToken() const
  {
    return this->Token;
  }
  const vtkFoamDict *GetUpperDictPtr() const
  {
    return this->UpperDictPtr;
  }

  // lookup by literal search
  vtkFoamEntry *LookupLiteral(const vtkStdString& keyword) const
  {
    if (this->Token.GetType() == vtkFoamToken::UNDEFINED)
      {
      for (size_t entryI = 0; entryI < this->Superclass::size(); entryI++)
        {
        // attempt literal match
        if (this->operator[](entryI)->GetKeyword() == keyword)
          {
          // found
          return this->operator[](entryI);
          }
        }
      }

    // not found
    return NULL;
  }

  // lookup by literal and regex search
  vtkFoamEntry *LookupRegEx(const vtkStdString& keyword) const
  {
    if (this->Token.GetType() == vtkFoamToken::UNDEFINED)
      {
      // first attempt a literal match regardless of the keywords being
      // matched containing regular expressions or not
      vtkFoamEntry *entry = this->LookupLiteral(keyword);
      if (entry == NULL)
        {
        // regex is matched in descending order
        for (int entryI = static_cast<int>(this->Superclass::size()) - 1;
            entryI >= 0; entryI--)
          {
          // attempt regex match
          if (this->operator[](entryI)->GetKeyword().RegExMatch(keyword))
            {
            // found
            entry = this->operator[](entryI);
            break;
            }
          }
        }
      return entry;
      }

    return NULL;
  }

  // return NULL if the found entry contains no value
  vtkFoamEntry *Lookup(const vtkStdString &keyword) const
  {
    vtkFoamEntry *entry = this->LookupRegEx(keyword);
    return (entry != NULL && entry->size() > 0) ? entry : NULL;
  }

  // reads a FoamFile or a subdictionary. if the stream to be read is
  // a subdictionary the preceding '{' is assumed to have already been
  // thrown away.
  bool Read(vtkFoamIOobject& io, const bool isSubDictionary = false,
      const vtkFoamToken& firstToken = vtkFoamToken())
  {
    try
      {
      bool isRegExKeyword = true;
      vtkFoamToken currToken;
      if (firstToken.GetType() == vtkFoamToken::UNDEFINED)
        {
        // read the first token
        if (!io.Read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }

        if (isSubDictionary)
          {
          // the following if clause is for an exceptional expression
          // of `LABEL{LABELorSCALAR}' without type prefix
          // (e. g. `2{-0}' in mixedRhoE B.C. in
          // rhopSonicFoam/shockTube)
          if (currToken.GetType() == vtkFoamToken::LABEL
              || currToken.GetType() == vtkFoamToken::SCALAR)
            {
            this->Token = currToken;
            io.ReadExpecting('}');
            return true;
            }
          // return as empty dictionary
          else if (currToken == '}')
            {
            return true;
            }
          }
        else
          {
          // list of dictionaries is read as a usual dictionary with regex
          // turned off: polyMesh/boundary, point/face/cell-Zones
          if (currToken.GetType() == vtkFoamToken::LABEL)
            {
            io.ReadExpecting('(');
            if (currToken.To<int>() > 0)
              {
              if (!io.Read(currToken))
                {
                throw vtkFoamError() << "Unexpected EOF";
                }
              // continue to read as a usual dictionary with regex turned off
              isRegExKeyword = false;
              }
            else // return as empty dictionary
              {
              io.ReadExpecting(')');
              return true;
              }
            }
          // some boundary files does not have the number of boundary
          // patches (e.g. settlingFoam/tank3D). In this case the file is
          // needed to be explicitly read as a dictionary.
          else if (currToken == '('
              && io.GetClassName() == "polyBoundaryMesh") // polyMesh/boundary
            {
            if (!io.Read(currToken)) // read the first keyword
              {
              throw vtkFoamError() << "Unexpected EOF";
              }
            if (currToken == ')') // return as empty dictionary
              {
              return true;
              }
            // continue to read as a usual dictionary with regex turned off
            isRegExKeyword = false;
            }
          }
        }
      // if firstToken is given as either a word or a string read the
      // following stream as subdictionary to constitute a list of
      // dictionaries. Both word and string are allowed as keyword,
      // but string does not work as regex.
      else if (firstToken.IsWordOrString())
        {
        this->Superclass::push_back(new vtkFoamEntry(this));
        this->Superclass::back()->SetKeyword(firstToken, isRegExKeyword);
        this->Superclass::back()->ReadDictionary(io);
        if (!io.Read(currToken) || currToken == '}' || currToken == ')')
          {
          return true;
          }
        }
      else // quite likely an identifier
        {
        currToken = firstToken;
        }

      if (currToken == ';' || currToken.IsWordOrString()
          || currToken.GetType() == vtkFoamToken::IDENTIFIER)
        {
        // general dictionary
        do
          {
          if (currToken.IsWordOrString())
            {
            vtkFoamEntry *previousEntry;
            if (isRegExKeyword && (previousEntry
                = this->LookupLiteral(currToken.ToStdString())) != NULL)
              {
              if (io.GetInputMode() == vtkFoamFile::INPUT_MODE_MERGE)
                {
                if (previousEntry->FirstValue().GetType()
                    == vtkFoamToken::DICTIONARY)
                  {
                  io.ReadExpecting('{');
                  previousEntry->FirstValue().Dictionary().Read(io, true);
                  }
                else
                  {
                  previousEntry->Clear();
                  previousEntry->Read(io);
                  }
                }
              else if (io.GetInputMode() == vtkFoamFile::INPUT_MODE_OVERWRITE)
                {
                previousEntry->Clear();
                previousEntry->Read(io);
                }
              else if (io.GetInputMode() == vtkFoamFile::INPUT_MODE_PROTECT)
                {
                // the contents of the entry is discarded
                vtkFoamEntry entry(this);
                entry.Read(io);
                }
              else // INPUT_MODE_ERROR || INPUT_MODE_WARN
                {
                // "#inputMode warn" just doesn't skip the duplicated entry but
                // in fact stops reading the remaining dictionary
                throw vtkFoamError() << "Found duplicated entries with keyword "
                  << currToken.ToStdString();
                }
              }
            else
              {
              this->Superclass::push_back(new vtkFoamEntry(this));
              this->Superclass::back()->SetKeyword(currToken, isRegExKeyword);
              this->Superclass::back()->Read(io);
              }

            // the "FoamFile" keyword of type word must be searched by exact
            // match (no need for regex search)
            if(currToken == "FoamFile")
              {
              // delete the FoamFile header subdictionary entry
              delete this->Superclass::back();
              this->Superclass::pop_back();
              }
            // the "include" keyword of type word or string must be searched by
            // exact match (no need for regex search)
            else if (currToken.ToStdString() == "include")
              {
              // include the named file. The name must be of type string.
              // Exiting the included file at EOF will be handled automatically
              // by vtkFoamFile::closeIncludedFile()
                if (this->Superclass::back()->size() == 0
                    || this->Superclass::back()->FirstValue().GetType()
                    != vtkFoamToken::STRING)
                {
                throw vtkFoamError()
                << "Expected string as the file name to be included, found "
                << this->Superclass::back()->size() == 0 ? vtkFoamToken()
                : this->Superclass::back()->FirstValue();
                }
              const vtkStdString includeFileName(
                  this->Superclass::back()->ToString());
              delete this->Superclass::back();
              this->Superclass::pop_back();
              io.IncludeFile(includeFileName, io.GetFilePath(), false);
              }
            }
          else if (currToken.GetType() == vtkFoamToken::IDENTIFIER)
            {
            // substitute identifier
            const vtkStdString identifier(currToken.ToStdString());

            for (const vtkFoamDict *uDictPtr = this; uDictPtr != NULL;)
              {
              const vtkFoamEntry *identifiedEntry
                  = uDictPtr->LookupLiteral(identifier); // do a regex search
              if (identifiedEntry != NULL)
                {
                if (identifiedEntry->size() == 0
                    || identifiedEntry->FirstValue().GetType()
                    != vtkFoamToken::DICTIONARY)
                  {
                  throw vtkFoamError()
                  << "Expected dictionary for substituting entry "
                  << identifier;
                  }
                const vtkFoamDict& identifiedDict
                = identifiedEntry->FirstValue().Dictionary();
                for (size_t entryI = 0; entryI < identifiedDict.size(); entryI++)
                  {
                  // I think #inputMode handling should be done here
                  // as well, but the genuine FoamFile parser for OF
                  // 1.5 does not seem to be doing it.
                  this->Superclass::push_back(
                      new vtkFoamEntry(*identifiedDict[entryI], this));
                  }
                break;
                }
              else
                {
                uDictPtr = uDictPtr->GetUpperDictPtr();
                if (uDictPtr == NULL)
                  {
                  // if no entry with the identifier as keyowrd is found
                  // interpret the identifier as keyword
                  this->Superclass::push_back(new vtkFoamEntry(this));
                  vtkFoamToken keywordToken;
                  keywordToken.SetWord("$" + identifier);
                  this->Superclass::back()->SetKeyword(keywordToken, isRegExKeyword);
                  this->Superclass::back()->Read(io);
                  }
                }
              }
            }
          // skip empty entry only with ';'
          }while (io.Read(currToken) && (currToken.IsWordOrString()
                || currToken.GetType() == vtkFoamToken::IDENTIFIER
                || currToken == ';'));

        if (currToken.GetType() == vtkFoamToken::TOKEN_ERROR || currToken == '}'
            || currToken == ')')
          {
          return true;
          }
        throw vtkFoamError()
        << "Expected keyword, closing brace, ';' or EOF, found " << currToken;
        }
      throw vtkFoamError() << "Expected keyword or identifier, found "
      << currToken;
      }
    catch (vtkFoamError& e)
      {
      if (isSubDictionary)
        {
        throw;
        }
      else
        {
        io.SetError(e);
        return false;
        }
      }
  }
};

void vtkFoamIOobject::ReadHeader()
{
  // the "FoamFile" keyword must be searched by exact match (no need
  // for regex search)
  this->Superclass::ReadExpecting("FoamFile");
  this->Superclass::ReadExpecting('{');

  vtkFoamDict headerDict;
  // throw exception in case of error
  headerDict.Read(*this, true, vtkFoamToken());

  const vtkFoamEntry *formatEntry = headerDict.Lookup("format");
  if (formatEntry == NULL || !formatEntry->FirstValue().IsWordOrString())
    {
    throw vtkFoamError()
    << "valid format entry (binary/ascii) not found in FoamFile header";
    }
  // case does matter (e. g. "BINARY" is treated as ascii)
  // cf. src/OpenFOAM/db/IOstreams/IOstreams/IOstream.C
  this->Format = (formatEntry->ToStdString() == "binary" ? BINARY : ASCII);

  const vtkFoamEntry *classEntry = headerDict.Lookup("class");
  if (classEntry == NULL || !classEntry->FirstValue().IsWordOrString())
    {
    throw vtkFoamError() << "valid class name not found in FoamFile header";
    }
  this->HeaderClassName = classEntry->ToStdString();

  const vtkFoamEntry *objectEntry = headerDict.Lookup("object");
  if (objectEntry == NULL || !objectEntry->FirstValue().IsWordOrString())
    {
    throw vtkFoamError() << "valid object name not found in FoamFile header";
    }
  this->ObjectName = objectEntry->ToStdString();
}

vtkFoamEntryValue::vtkFoamEntryValue(
    vtkFoamEntryValue& value, const vtkFoamEntry *upperEntryPtr) :
  vtkFoamToken(value), IsUniform(value.GetIsUniform()), Managed(true),
      UpperEntryPtr(upperEntryPtr)
{
  switch (this->Superclass::Type)
    {
    case VECTORLIST:
#if VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS
      {
      vtkFloatArray *fa = vtkFloatArray::SafeDownCast(value.ToVTKObject());
      if(fa->GetNumberOfComponents() == 6)
        {
        // create deepcopies for vtkObjects to avoid duplicated
        // mainpulation of symmTensor components
        vtkFloatArray *newfa = vtkFloatArray::New();
        newfa->DeepCopy(fa);
        this->Superclass::VtkObjectPtr = newfa;
        break;
        }
      }
#endif
    case LABELLIST:
    case SCALARLIST:
    case STRINGLIST:
      this->Superclass::VtkObjectPtr = value.ToVTKObject();
      this->Superclass::VtkObjectPtr->Register(0);
      break;
    case LABELLISTLIST:
      this->LabelListListPtr = new vtkFoamIntVectorVector(*value.LabelListListPtr);
      break;
    case ENTRYVALUELIST:
      {
      const size_t nValues = value.EntryValuePtrs->size();
      this->EntryValuePtrs = new vtkstd::vector<vtkFoamEntryValue*>(nValues);
      for (size_t valueI = 0; valueI < nValues; valueI++)
        {
        this->EntryValuePtrs->operator[](valueI) = new vtkFoamEntryValue(
            *value.EntryValuePtrs->operator[](valueI), upperEntryPtr);
        }
      }
      break;
    case DICTIONARY:
      // UpperEntryPtr is null when called from vtkFoamDict constructor
      if (this->UpperEntryPtr != NULL)
        {
        this->DictPtr = new vtkFoamDict(*value.DictPtr,
            this->UpperEntryPtr->GetUpperDictPtr());
        }
      else
        {
        this->DictPtr = NULL;
        }
      break;
    case EMPTYLIST:
      break;
      // required to suppress the 'enumeration value not handled' warning by
      // g++ when compiled with -Wall
    default:
      break;
    }
}

void vtkFoamEntryValue::Clear()
{
  this->IsUniform = UNDEFINED;
  if (this->Managed)
    {
    switch (this->Superclass::Type)
      {
      case LABELLIST:
      case SCALARLIST:
      case VECTORLIST:
      case STRINGLIST:
        this->VtkObjectPtr->Delete();
        break;
      case LABELLISTLIST:
        delete this->LabelListListPtr;
        break;
      case ENTRYVALUELIST:
        for (size_t valueI = 0; valueI < this->EntryValuePtrs->size() ; valueI++)
          {
          delete this->EntryValuePtrs->operator[](valueI);
          }
        delete this->EntryValuePtrs;
        break;
      case DICTIONARY:
        delete this->DictPtr;
        break;
        // required to suppress the 'enumeration value not handled' warning by
        // g++ when compiled with -Wall
      default:
        break;
      }
    }
}

// general-purpose list reader - guess the type of the list and read
// it. only supports ascii format and assumes the preceding '(' has
// already been thrown away.  the reader supports nested list with
// variable lengths (e. g. `((token token) (token token token)).'
// also supports compound of tokens and lists (e. g. `((token token)
// token)') only if a list comes as the first value.
void vtkFoamEntryValue::ReadList(vtkFoamIOobject& io)
{
  vtkFoamToken currToken;
  io.Read(currToken);

  // initial guess of the list type
  if (currToken.GetType() == this->Superclass::LABEL)
    {
    // if the first token is of type LABEL it might be either an element of
    // a labelList or the size of a sublist so proceed to the next token
    vtkFoamToken nextToken;
    if (!io.Read(nextToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if (nextToken.GetType() == this->Superclass::LABEL)
      {
      this->Superclass::LabelListPtr = vtkIntArray::New();
      this->Superclass::LabelListPtr->InsertNextValue(currToken.To<int>());
      this->Superclass::LabelListPtr->InsertNextValue(nextToken.To<int>());
      this->Superclass::Type = LABELLIST;
      }
    else if (nextToken.GetType() == this->Superclass::SCALAR)
      {
      this->Superclass::ScalarListPtr = vtkFloatArray::New();
      this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
      this->Superclass::ScalarListPtr->InsertNextValue(nextToken.To<float>());
      this->Superclass::Type = SCALARLIST;
      }
    else if (nextToken == '(') // list of list: read recursively
      {
      this->Superclass::EntryValuePtrs = new vtkstd::vector<vtkFoamEntryValue*>;
      this->Superclass::EntryValuePtrs->push_back(new vtkFoamEntryValue(
          this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      this->Superclass::Type = ENTRYVALUELIST;
      }
    else if (nextToken == ')') // list with only one label element
      {
      this->Superclass::LabelListPtr = vtkIntArray::New();
      this->Superclass::LabelListPtr->SetNumberOfValues(1);
      this->Superclass::LabelListPtr->SetValue(0, currToken.To<int>());
      this->Superclass::Type = LABELLIST;
      return;
      }
    else
      {
      throw vtkFoamError() << "Expected number, '(' or ')', found "
      << nextToken;
      }
    }
  else if (currToken.GetType() == this->Superclass::SCALAR)
    {
    this->Superclass::ScalarListPtr = vtkFloatArray::New();
    this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
    this->Superclass::Type = SCALARLIST;
    }
  // if the first word is a string we have to read another token to determine
  // if the first word is a keyword for the following dictionary
  else if (currToken.IsWordOrString())
    {
    vtkFoamToken nextToken;
    if (!io.Read(nextToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if (nextToken.IsWordOrString()) // list of strings
      {
      this->Superclass::StringListPtr = vtkStringArray::New();
      this->Superclass::StringListPtr->InsertNextValue(currToken.ToStdString());
      this->Superclass::StringListPtr->InsertNextValue(nextToken.ToStdString());
      this->Superclass::Type = STRINGLIST;
      }
    // dictionary with the already read word or string as the first keyword
    else if (nextToken == '{')
      {
      this->ReadDictionary(io, currToken);
      // the dictionary read as list has the entry terminator ';' so
      // we have to skip it
      return;
      }
    else if (nextToken == ')') // list with only one string element
      {
      this->Superclass::StringListPtr = vtkStringArray::New();
      this->Superclass::StringListPtr->SetNumberOfValues(1);
      this->Superclass::StringListPtr->SetValue(0, currToken.ToStdString());
      this->Superclass::Type = STRINGLIST;
      return;
      }
    else
      {
      throw vtkFoamError() << "Expected string, '{' or ')', found "
      << nextToken;
      }
    }
  // list of lists or dictionaries: read recursively
  else if (currToken == '(' || currToken == '{')
    {
    this->Superclass::EntryValuePtrs = new vtkstd::vector<vtkFoamEntryValue*>;
    this->Superclass::EntryValuePtrs->push_back(new vtkFoamEntryValue(
        this->UpperEntryPtr));
    if(currToken == '(')
      {
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      }
    else // currToken == '{'
      {
      this->Superclass::EntryValuePtrs->back()->ReadDictionary(io, vtkFoamToken());
      }
    // read all the following values as arbitrary entryValues
    // the alphaContactAngle b.c. in multiphaseInterFoam/damBreak4phase
    // reaquires this treatment (reading by readList() is not enough)
    do
      {
      this->Superclass::EntryValuePtrs->push_back(new vtkFoamEntryValue(
          this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->Read(io);
      } while (*this->Superclass::EntryValuePtrs->back() != ')'
        && *this->Superclass::EntryValuePtrs->back() != '}'
        && *this->Superclass::EntryValuePtrs->back() != ';');

    if (*this->Superclass::EntryValuePtrs->back() != ')')
      {
      throw vtkFoamError() << "Expected ')' before "
          << *this->Superclass::EntryValuePtrs->back();
      }

    // delete ')'
    delete this->Superclass::EntryValuePtrs->back();
    this->EntryValuePtrs->pop_back();
    this->Superclass::Type = ENTRYVALUELIST;
    return;
    }
  else if (currToken == ')') // empty list
    {
    this->Superclass::Type = EMPTYLIST;
    return;
    }
  // FIXME: may (or may not) need identifier handling

  while (io.Read(currToken) && currToken != ')')
    {
    if (this->Superclass::Type == LABELLIST)
      {
      if (currToken.GetType() == this->Superclass::SCALAR)
        {
        // switch to scalarList
        // LabelListPtr and ScalarListPtr are packed into a single union so
        // we need a temprary pointer
        vtkFloatArray* slPtr = vtkFloatArray::New();
        const int size = this->Superclass::LabelListPtr->GetNumberOfTuples();
        slPtr->SetNumberOfValues(size + 1);
        for (int i = 0; i < size; i++)
          {
          slPtr->SetValue(i,
              static_cast<float>(this->Superclass::LabelListPtr->GetValue(i)));
          }
        this->LabelListPtr->Delete();
        slPtr->SetValue(size, currToken.To<float>());
        // copy after LabelListPtr is deleted
        this->Superclass::ScalarListPtr = slPtr;
        this->Superclass::Type = SCALARLIST;
        }
      else if (currToken.GetType() == this->Superclass::LABEL)
        {
        this->Superclass::LabelListPtr->InsertNextValue(currToken.To<int>());
        }
      else
        {
        throw vtkFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::SCALARLIST)
      {
      if (currToken.Is<float>())
        {
        this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
        }
      else
        {
        throw vtkFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::STRINGLIST)
      {
      if (currToken.IsWordOrString())
        {
        this->Superclass::StringListPtr->InsertNextValue(currToken.ToStdString());
        }
      else
        {
        throw vtkFoamError() << "Expected a string, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::ENTRYVALUELIST)
      {
      if (currToken.GetType() == this->Superclass::LABEL)
        {
        // skip the number of elements to make things simple
        if (!io.Read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }
        }
      if (currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      this->Superclass::EntryValuePtrs->push_back(new vtkFoamEntryValue(this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      }
    else
      {
      throw vtkFoamError() << "Unexpected token " << currToken;
      }
    }

  if (this->Superclass::Type == this->Superclass::LABELLIST)
    {
    this->Superclass::LabelListPtr->Squeeze();
    }
  else if (this->Superclass::Type == this->Superclass::SCALARLIST)
    {
    this->Superclass::ScalarListPtr->Squeeze();
    }
  else if (this->Superclass::Type == this->Superclass::STRINGLIST)
    {
    this->Superclass::StringListPtr->Squeeze();
    }
}

// a list of dictionaries is actually read as a dictionary
void vtkFoamEntryValue::ReadDictionary(vtkFoamIOobject& io,
    const vtkFoamToken& firstKeyword)
{
  this->Superclass::DictPtr = new vtkFoamDict(this->UpperEntryPtr->GetUpperDictPtr());
  this->Superclass::Type = this->Superclass::DICTIONARY;
  this->Superclass::DictPtr->Read(io, true, firstKeyword);
}

// guess the type of the given entry value and read it
void vtkFoamEntryValue::Read(vtkFoamIOobject& io)
{
  vtkFoamToken currToken;
  if (!io.Read(currToken))
    {
    throw vtkFoamError() << "Unexpected EOF";
    }

  if (currToken == '{')
    {
    this->ReadDictionary(io, vtkFoamToken());
    return;
    }
  // for reading sublist from vtkFoamEntryValue::readList() or there
  // are cases where lists without the (non)uniform keyword appear
  // (e. g. coodles/pitsDaily/0/U, uniformFixedValue b.c.)
  else if (currToken == '(')
    {
    this->ReadList(io);
    return;
    }
  else if (currToken == '[')
    {
    this->ReadDimensionSet(io);
    return;
    }
  else if (currToken == "uniform")
    {
    if (!io.Read(currToken))
      {
      throw vtkFoamError()
      << "Expected a uniform value or a list, found unexpected EOF";
      }
    if (currToken == '(')
      {
      this->ReadList(io);
      }
    else if (currToken.GetType() == this->Superclass::LABEL
        || currToken.GetType() == this->Superclass::SCALAR
        || currToken.IsWordOrString()
        || currToken.GetType() == this->Superclass::IDENTIFIER)
      {
      this->Superclass::operator=(currToken);
      }
    else // unexpected punctuation token
      {
      throw vtkFoamError() << "Expected number, string or (, found "
      << currToken;
      }
    this->IsUniform = UNIFORM;
    }
  else if (currToken == "nonuniform")
    {
    if (!io.Read(currToken))
      {
      throw vtkFoamError() << "Expected list type specifier, found EOF";
      }
    this->IsUniform = NONUNIFORM;
    if (currToken == "List<scalar>")
      {
      this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(io);
      }
    else if (currToken == "List<sphericalTensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 1, false> >(io);
      }
    else if (currToken == "List<vector>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 3, false> >(io);
      }
    else if (currToken == "List<symmTensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 6, false> >(io);
      }
    else if (currToken == "List<tensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 9, false> >(io);
      }
    // List<bool> is read as List<label>
    else if (currToken == "List<label>" || currToken == "List<bool>")
      {
      this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
      }
    // an empty list doesn't have a list type specifier
    else if (currToken.GetType() == this->Superclass::LABEL
        && currToken.To<int>() == 0)
      {
      this->Superclass::Type = this->Superclass::EMPTYLIST;
      if(io.GetFormat() == vtkFoamIOobject::ASCII)
        {
        io.ReadExpecting('(');
        io.ReadExpecting(')');
        }
      }
    else if (currToken.GetType() == this->Superclass::IDENTIFIER)
      {
      this->Superclass::operator=(currToken);
      }
    else
      {
      throw vtkFoamError() << "Unsupported nonuniform list type " << currToken;
      }
    }
  // the followings read nonuniform lists without setting isUniform to
  // NONUNIFORM
  else if (currToken == "List<scalar>")
    {
    this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(io);
    }
  else if (currToken == "List<sphericalTensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 1, false> >(io);
    }
  else if (currToken == "List<vector>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 3, false> >(io);
    }
  else if (currToken == "List<symmTensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 6, false> >(io);
    }
  else if (currToken == "List<tensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 9, false> >(io);
    }
  // zones have list without a uniform/nonuniform keyword
  // List<bool> is read as List<label>
  // (e. g. flipMap entry in faceZones)
  else if (currToken == "List<label>" || currToken == "List<bool>")
    {
    this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
    }
  else if (currToken.GetType() == this->Superclass::PUNCTUATION
      || currToken.GetType() == this->Superclass::LABEL || currToken.GetType()
      == this->Superclass::SCALAR || currToken.IsWordOrString()
      || currToken.GetType() == this->Superclass::IDENTIFIER)
    {
    this->Superclass::operator=(currToken);
    }
}

// read values of an entry
void vtkFoamEntry::Read(vtkFoamIOobject& io)
{
  for (;;)
    {
    this->Superclass::push_back(new vtkFoamEntryValue(this));
    this->Superclass::back()->Read(io);

    if (this->Superclass::size() >= 2)
      {
      vtkFoamEntryValue& secondLastValue =
          *this->Superclass::operator[](this->Superclass::size() - 2);
      if (secondLastValue.GetType() == vtkFoamToken::LABEL)
        {
        vtkFoamEntryValue& lastValue = *this->Superclass::back();

        // a zero-sized nonuniform list without prefixing "nonuniform"
        // keyword nor list type specifier (i. e. `0()';
        // e. g. simpleEngine/0/polyMesh/pointZones) requires special
        // care (one with nonuniform prefix is treated within
        // vtkFoamEntryValue::read()). still this causes errornous
        // behavior for `0 nonuniform 0()' but this should be extremely
        // rare
        if (lastValue.GetType() == vtkFoamToken::EMPTYLIST && secondLastValue
            == 0)
          {
          delete this->Superclass::back();
          this->Superclass::pop_back(); // delete the last value
          // mark new last value as empty
          this->Superclass::back()->SetEmptyList();
          }
        // for an exceptional expression of `LABEL{LABELorSCALAR}' without
        // type prefix (e. g. `2{-0}' in mixedRhoE B.C. in
        // rhopSonicFoam/shockTube)
        else if (lastValue.GetType() == vtkFoamToken::DICTIONARY)
          {
          if (lastValue.Dictionary().GetType() == vtkFoamToken::LABEL)
            {
            const int asize = secondLastValue.To<int>();
            // not using templated To<int>() for workarounding an error
            // on SunOS-CC5.6-dbg
            const int value = lastValue.Dictionary().GetToken().ToInt();
            // delete last two values
            delete this->Superclass::back();
            this->Superclass::pop_back();
            delete this->Superclass::back();
            this->Superclass::pop_back();
            // make new labelList
            this->Superclass::push_back(new vtkFoamEntryValue(this));
            this->Superclass::back()->MakeLabelList(value, asize);
            }
          else if (lastValue.Dictionary().GetType() == vtkFoamToken::SCALAR)
            {
            const int asize = secondLastValue.To<int>();
            // not using templated To<float>() for workarounding an error
            // on SunOS-CC5.6-dbg
            const float value = lastValue.Dictionary().GetToken().ToFloat();
            // delete last two values
            delete this->Superclass::back();
            this->Superclass::pop_back();
            delete this->Superclass::back();
            this->Superclass::pop_back();
            // make new labelList
            this->Superclass::push_back(new vtkFoamEntryValue(this));
            this->Superclass::back()->MakeScalarList(value, asize);
            }
          }
        }
      }

    if (this->Superclass::back()->GetType() == vtkFoamToken::IDENTIFIER)
      {
      // substitute identifier
      const vtkStdString identifier(this->Superclass::back()->ToStdString());
      vtkFoamEntryValue::uniformTypes isUniform
          = this->Superclass::back()->GetIsUniform();
      delete this->Superclass::back();
      this->Superclass::pop_back();

      for (const vtkFoamDict *uDictPtr = this->UpperDictPtr; uDictPtr != NULL;)
        {
        const vtkFoamEntry *identifiedEntry = uDictPtr->LookupLiteral(identifier);

        if (identifiedEntry != NULL)
          {
          if (isUniform != vtkFoamEntryValue::UNDEFINED
              && identifiedEntry->size() > 0
              && identifiedEntry->operator[](0)->GetIsUniform()
              != vtkFoamEntryValue::UNDEFINED)
            {
            throw vtkFoamError() << "duplicated list type specifiers (uniform/"
              "nonuniform) in the first substituted entry value of identifier $"
              << identifier;
            }

          for (size_t valueI = 0; valueI < identifiedEntry->size(); valueI++)
            {
            if (identifiedEntry->operator[](valueI)->GetType()
                == vtkFoamToken::DICTIONARY)
              {
              // in fact there's no problem with substituting
              // dictionary in this parser. Only for compatibility
              // with the genuine parser
              throw vtkFoamError() << "dictionary cannot be substituted as an "
                  "entry value. Try {$" << identifier << ";} instead of $"
                  << identifier << ".";
              }
            this->Superclass::push_back(new vtkFoamEntryValue(
                *identifiedEntry->operator[](valueI), this));
            }

          if (isUniform != vtkFoamEntryValue::UNDEFINED)
            {
            if (identifiedEntry->size() > 0)
              {
              this->Superclass::operator[](this->Superclass::size()
                  - identifiedEntry->size())->SetIsUniform(isUniform);
              }
            else
              {
              this->Superclass::push_back(new vtkFoamEntryValue(this));
              this->Superclass::back()->SetIsUniform(isUniform);
              }
            }
          break;
          }
        else
          {
          uDictPtr = uDictPtr->GetUpperDictPtr();
          if (uDictPtr == NULL)
            {
            // if the substituting entry is not found append the identifier as
            // word
            this->Superclass::push_back(new vtkFoamEntryValue(this));
            this->Superclass::back()->SetWord("$" + identifier);
            this->Superclass::back()->SetIsUniform(isUniform);
            }
          }
        }
      }
    else if (*this->Superclass::back() == ';')
      {
      delete this->Superclass::back();
      this->Superclass::pop_back();
      break;
      }
    else if (this->Superclass::back()->GetType() == vtkFoamToken::DICTIONARY)
      {
      // subdictionary is not suffixed by an entry terminator ';'
      break;
      }
    else if (*this->Superclass::back() == '}' || *this->Superclass::back()
        == ')')
      {
      throw vtkFoamError() << "Unmatched " << *this->Superclass::back();
      }
    }
}

//-----------------------------------------------------------------------------
// vtkOFFReaderPrivate static members
const char *vtkOFFReaderPrivate::InternalMeshIdentifier
    = "[InternalMesh]";
const char *vtkOFFReaderPrivate::SurfaceMeshIdentifier
    = "[InternalSurfaceMesh]";
const char *vtkOFFReaderPrivate::ProcessorPatchesIdentifier
    = "[ProcessorPatches]";

//-----------------------------------------------------------------------------
int vtkOFFReaderPrivate::GetBoundaryTypeProcessor()
{
  return static_cast<int>(vtkFoamBoundaryEntry::PROCESSOR);
}

//-----------------------------------------------------------------------------
int vtkOFFReaderPrivate::GetBoundaryTypeInternal()
{
  return static_cast<int>(vtkFoamBoundaryEntry::INTERNAL);
}

//-----------------------------------------------------------------------------
// vtkOFFReaderPrivate constructor and destructor
vtkOFFReaderPrivate::vtkOFFReaderPrivate()
{
  // DATA TIMES
  this->TimeStep = 0;
  this->TimeStepOld = -1;
  this->TimeValues = vtkDoubleArray::New();
  this->TimeNames = vtkStringArray::New();

  // selection
  this->InternalMeshSelectionStatus = 0;
  this->InternalMeshSelectionStatusOld = 0;

  // DATA COUNTS
  this->NumCells = 0;
  this->NumPoints = 0;
  this->NumInternalFaces = 0;

  this->VolFieldFiles = vtkStringArray::New();
  this->SurfaceFieldFiles = vtkStringArray::New();
  this->PointFieldFiles = vtkStringArray::New();
  this->LagrangianFieldFiles = vtkStringArray::New();
  this->PolyMeshPointsDir = vtkStringArray::New();
  this->PolyMeshFacesDir = vtkStringArray::New();

  // for creating cell-to-point translated data
  this->BoundaryPointMap = NULL;
  this->AllBoundaries = NULL;
  this->AllBoundariesPointMap = NULL;
  this->InternalPoints = NULL;

  // for caching mesh
  this->InternalMesh = NULL;
  this->SurfaceMesh = NULL;
  this->BoundaryMesh = NULL;
  this->BoundaryPointMap = NULL;
#if 0
  this->ReciprocalDelta = NULL;
#endif
  this->FaceOwner = NULL;
  this->ProcessorFaces = NULL;
  this->LagrangianMesh = NULL;
  this->PointZoneMesh = NULL;
  this->FaceZoneMesh = NULL;
  this->CellZoneMesh = NULL;

  // for decomposing polyhedra
  this->NumAdditionalCells = 0;
  this->AdditionalCellIds = NULL;
  this->NumAdditionalCells = NULL;
  this->AdditionalCellPoints = NULL;
}

//-----------------------------------------------------------------------------
vtkOFFReaderPrivate::~vtkOFFReaderPrivate()
{
  this->TimeValues->Delete();
  this->TimeNames->Delete();

  this->VolFieldFiles->Delete();
  this->SurfaceFieldFiles->Delete();
  this->PointFieldFiles->Delete();
  this->LagrangianFieldFiles->Delete();
  this->PolyMeshPointsDir->Delete();
  this->PolyMeshFacesDir->Delete();

  this->ClearMeshes();
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Case Path: "
      << (this->CasePath.length() ? this->CasePath.c_str() : "(none)") << endl;
  os << indent << "Region Name: "
      << (this->RegionName.length() ? this->RegionName.c_str() : "(none)")
      << endl;
  os << indent << "Processor Name: "
      << (this->ProcessorName.length() ? this->ProcessorName.c_str() : "(none)")
      << endl;
  if (this->TimeValues)
    {
    os << indent << "Number of Time Steps: "
        << this->TimeValues->GetNumberOfTuples() << endl;
    }
  os << indent << "Time Step: " << this->TimeStep << endl;
  os << indent << "Number of Cells: " << this->NumCells << endl;
  os << indent << "Number of Points: " << this->NumPoints << endl;
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::ClearInternalMeshes()
{
  if (this->FaceOwner != NULL)
    {
    this->FaceOwner->Delete();
    this->FaceOwner = NULL;
    }
  delete this->ProcessorFaces;
  this->ProcessorFaces = NULL;
  if (this->InternalMesh != NULL)
    {
    this->InternalMesh->Delete();
    this->InternalMesh = NULL;
    }
  if (this->SurfaceMesh != NULL)
    {
    this->SurfaceMesh->Delete();
    this->SurfaceMesh = NULL;
    }
  if (this->AdditionalCellIds != NULL)
    {
    this->AdditionalCellIds->Delete();
    this->AdditionalCellIds = NULL;
    }
  if (this->NumAdditionalCells != NULL)
    {
    this->NumAdditionalCells->Delete();
    this->NumAdditionalCells = NULL;
    }
  delete this->AdditionalCellPoints;
  this->AdditionalCellPoints = NULL;

  if (this->PointZoneMesh != NULL)
    {
    this->PointZoneMesh->Delete();
    this->PointZoneMesh = NULL;
    }
  if (this->FaceZoneMesh != NULL)
    {
    this->FaceZoneMesh->Delete();
    this->FaceZoneMesh = NULL;
    }
  if (this->CellZoneMesh != NULL)
    {
    this->CellZoneMesh->Delete();
    this->CellZoneMesh = NULL;
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::ClearBoundaryMeshes()
{
  if (this->BoundaryMesh != NULL)
    {
    this->BoundaryMesh->Delete();
    this->BoundaryMesh = NULL;
    }

  delete this->BoundaryPointMap;
  this->BoundaryPointMap = NULL;
#if 0
  delete this->ReciprocalDelta;
  this->ReciprocalDelta = NULL;
#endif

  if (this->InternalPoints != NULL)
    {
    this->InternalPoints->Delete();
    this->InternalPoints = NULL;
    }
  if (this->AllBoundaries != NULL)
    {
    this->AllBoundaries->Delete();
    this->AllBoundaries = NULL;
    }
  if (this->AllBoundariesPointMap != NULL)
    {
    this->AllBoundariesPointMap->Delete();
    this->AllBoundariesPointMap = NULL;
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::ClearLagrangianMeshes()
{
  if (this->LagrangianMesh)
    {
    this->LagrangianMesh->Delete();
    this->LagrangianMesh = 0;
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::ClearMeshes()
{
  this->ClearInternalMeshes();
  this->ClearBoundaryMeshes();
  this->ClearLagrangianMeshes();
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::SetTimeValue(const double requestedTime)
{
  const int nTimeValues = this->TimeValues->GetNumberOfTuples();
  if (nTimeValues > 0)
    {
    int minTimeI = 0;
    double minTimeDiff = fabs(this->TimeValues->GetValue(0) - requestedTime);
    for (int timeI = 1; timeI < nTimeValues; timeI++)
      {
      const double timeDiff(fabs(this->TimeValues->GetValue(timeI)
          - requestedTime));
      if (timeDiff < minTimeDiff)
        {
        minTimeI = timeI;
        minTimeDiff = timeDiff;
        }
      }
    this->SetTimeStep(minTimeI); // set Modified() if TimeStep changed
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::SetupInformation(const vtkStdString &casePath,
    const vtkStdString &regionName, const vtkStdString &procName,
    vtkOFFReaderPrivate *master)
{
  // copy parent, path and timestep information from master
  this->CasePath = casePath;
  this->RegionName = regionName;
  this->ProcessorName = procName;
  this->Parent = master->Parent;
  this->TimeValues->Delete();
  this->TimeValues = master->TimeValues;
  this->TimeValues->Register(0);
  this->TimeNames->Delete();
  this->TimeNames = master->TimeNames;
  this->TimeNames->Register(0);

  this->PopulatePolyMeshDirArrays();
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::GetFieldNames(const vtkStdString &tempPath,
    const bool isLagrangian, vtkStringArray *cellObjectNames,
    vtkStringArray *surfaceObjectNames, vtkStringArray *pointObjectNames)
{
  // open the directory and get num of files
  vtkDirectory *directory = vtkDirectory::New();
  if (!directory->Open(tempPath.c_str()))
    {
    // no data
    directory->Delete();
    return;
    }

  // loop over all files and locate valid fields
  int nFieldFiles = directory->GetNumberOfFiles();
  for (int j = 0; j < nFieldFiles; j++)
    {
    const vtkStdString fieldFile(directory->GetFile(j));
    const size_t len = fieldFile.length();

    // excluded extensions cf. src/OpenFOAM/OSspecific/Unix/Unix.C
    if (!directory->FileIsDirectory(fieldFile.c_str()) && fieldFile.substr(len
        - 1) != "~" && (len < 4 || (fieldFile.substr(len - 4) != ".bak"
        && fieldFile.substr(len - 4) != ".BAK" && fieldFile.substr(len - 4)
        != ".old")) && (len < 5 || fieldFile.substr(len - 5) != ".save"))
      {
      vtkFoamIOobject io(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      if (io.Open(tempPath + "/" + fieldFile)) // file exists and readable
        {
        const vtkStdString& cn = io.GetClassName();
        if (isLagrangian)
          {
          if (cn == "labelField" || cn == "scalarField" || cn == "vectorField"
              || cn == "sphericalTensorField" || cn == "symmTensorField" || cn
              == "tensorField")
            {
            // real file name
            this->LagrangianFieldFiles->InsertNextValue(fieldFile);
            // object name
            pointObjectNames->InsertNextValue(io.GetObjectName());
            }
          }
        else
          {
          if (cn == "volScalarField" || cn == "surfaceScalarField"
              || cn == "pointScalarField" || cn == "volVectorField"
              || cn == "surfaceVectorField" || cn == "pointVectorField"
              || cn == "volSphericalTensorField"
              || cn == "surfaceSphericalTensorField"
              || cn == "pointSphericalTensorField"
              || cn == "volSymmTensorField" || cn == "surfaceSymmTensorField"
              || cn == "pointSymmTensorField" || cn == "volTensorField"
              || cn == "surfaceTensorField" || cn == "pointTensorField"
              || cn == "volScalarField::DimensionedInternalField"
              || cn == "volVectorField::DimensionedInternalField"
              || cn == "volSphericalTensorField::DimensionedInternalField"
              || cn == "volSymmTensorField::DimensionedInternalField"
              || cn == "volTensorField::DimensionedInternalField")
            {
            if (cn.substr(0, 3) == "vol")
              {
              // real file name
              this->VolFieldFiles->InsertNextValue(fieldFile);
              // object name
              cellObjectNames->InsertNextValue(io.GetObjectName());
              }
            else if (cn.substr(0, 7) == "surface")
              {
              // real file name
              this->SurfaceFieldFiles->InsertNextValue(fieldFile);
              // object name
              surfaceObjectNames->InsertNextValue(io.GetObjectName());
              }
            else
              {
              this->PointFieldFiles->InsertNextValue(fieldFile);
              pointObjectNames->InsertNextValue(io.GetObjectName());
              }
            }
          }
        io.Close();
        }
      }
    }
  // inserted objects are squeezed later in SortFieldFiles()
  directory->Delete();
}

//-----------------------------------------------------------------------------
// locate laglangian clouds
void vtkOFFReaderPrivate::LocateLagrangianClouds(
    vtkStringArray *lagrangianObjectNames, const vtkStdString &timePath)
{
  vtkDirectory *directory = vtkDirectory::New();
  if (directory->Open((timePath + this->RegionPath() + "/lagrangian").c_str()))
    {
    // search for sub-clouds (OF 1.5 format)
    const int nFiles = directory->GetNumberOfFiles();
    bool isSubCloud = false;
    for (int fileI = 0; fileI < nFiles; fileI++)
      {
      const vtkStdString fileNameI(directory->GetFile(fileI));
      if (fileNameI != "." && fileNameI != ".."
          && directory->FileIsDirectory(fileNameI.c_str()))
        {
        vtkFoamIOobject io(this->CasePath,
            this->Parent->GetIsSinglePrecisionBinary() != 0);
        const vtkStdString subCloudName(this->RegionPrefix() + "lagrangian/"
            + fileNameI);
        const vtkStdString subCloudFullPath(timePath + "/" + subCloudName);
        // lagrangian positions. there are many concrete class names
        // e. g. Cloud<parcel>, basicKinematicCloud etc.
        if ((io.Open(subCloudFullPath + "/positions")
            || io.Open(subCloudFullPath + "/positions.gz")) && io.GetClassName().find("Cloud") != vtkStdString::npos && io.GetObjectName()
            == "positions")
          {
          isSubCloud = true;
          // a lagrangianPath has to be in a bit different format from
          // subCloudName to make the "lagrangian" reserved path
          // component and a mesh region with the same name
          // distinguishable later
          const vtkStdString subCloudPath(this->RegionName + "/lagrangian/"
              + fileNameI);
          if (this->Parent->LagrangianPaths->LookupValue(subCloudPath) == -1)
            {
            this->Parent->LagrangianPaths->InsertNextValue(subCloudPath);
            }
          this->GetFieldNames(subCloudFullPath, true, NULL, NULL,
              lagrangianObjectNames);
          this->Parent->PatchDataArraySelection->AddArray(subCloudName.c_str());
          }
        }
      }
    // if there's no sub-cloud then OF < 1.5 format
    if (!isSubCloud)
      {
      vtkFoamIOobject io(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      const vtkStdString cloudName(this->RegionPrefix() + "lagrangian");
      const vtkStdString cloudFullPath(timePath + "/" + cloudName);
      if ((io.Open(cloudFullPath + "/positions") || io.Open(cloudFullPath
          + "/positions.gz")) && io.GetClassName().find("Cloud") != vtkStdString::npos && io.GetObjectName()
          == "positions")
        {
        const vtkStdString cloudPath(this->RegionName + "/lagrangian");
        if (this->Parent->LagrangianPaths->LookupValue(cloudPath) == -1)
          {
          this->Parent->LagrangianPaths->InsertNextValue(cloudPath);
          }
        this->GetFieldNames(cloudFullPath, true, NULL, NULL,
            lagrangianObjectNames);
        this->Parent->PatchDataArraySelection->AddArray(cloudName.c_str());
        }
      }
    this->Parent->LagrangianPaths->Squeeze();
    }
  directory->Delete();
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::SortFieldFiles(vtkStringArray *selections,
    vtkStringArray *files, vtkStringArray *objects)
{
  objects->Squeeze();
  files->Squeeze();
  vtkSortDataArray::Sort(objects, files);
  for (int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->InsertNextValue(objects->GetValue(nameI));
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
// create field data lists and cell/point array selection lists
int vtkOFFReaderPrivate::MakeMetaDataAtTimeStep(
    vtkStringArray *cellSelectionNames, vtkStringArray *surfaceSelectionNames,
    vtkStringArray *pointSelectionNames,
    vtkStringArray *lagrangianSelectionNames, const bool listNextTimeStep)
{
  // Read the patches from the boundary file into selection array
  if (this->PolyMeshFacesDir->GetValue(this->TimeStep)
      != this->BoundaryDict.TimeDir
      || this->Parent->PatchDataArraySelection->GetMTime()
          != this->Parent->PatchSelectionMTimeOld)
    {
    this->BoundaryDict.clear();
    this->BoundaryDict.TimeDir
        = this->PolyMeshFacesDir->GetValue(this->TimeStep);
    this->BoundaryDict.UpperProcRanges.clear();
    this->BoundaryDict.LowestUpperProcFaceNo = VTK_INT_MAX;
    this->BoundaryDict.ProcBoundaries.clear();

    vtkFoamDict *boundaryDict = this->GatherBlocks("boundary", this->TimeStep);
    if (boundaryDict == NULL && listNextTimeStep
        && this->TimeValues->GetNumberOfTuples() >= 2 && this->TimeStep == 0)
      {
      this->BoundaryDict.TimeDir
          = this->PolyMeshFacesDir->GetValue(1);
      boundaryDict = this->GatherBlocks("boundary", 1);
      }
    if (boundaryDict != NULL)
      {
      // Add the internal mesh by default always
      const vtkStdString internalMeshName(this->RegionPrefix()
          + vtkOFFReaderPrivate::InternalMeshIdentifier);
      if (this->RegionName == ""
          && this->Parent->Readers->GetNumberOfItems() > 1
          && !this->Parent->PatchDataArraySelection
          ->ArrayExists(internalMeshName.c_str()))
        {
        this->Parent->PatchDataArraySelection->DisableArray(
            internalMeshName.c_str());
        }
      else
        {
        this->Parent->PatchDataArraySelection
            ->AddArray(internalMeshName.c_str());
        }
      this->InternalMeshSelectionStatus
          = this->Parent->GetPatchArrayStatus(internalMeshName.c_str());

      // Add the surface mesh by default always
      const vtkStdString surfaceMeshName(this->RegionPrefix()
          + vtkOFFReaderPrivate::SurfaceMeshIdentifier);
      if (!this->Parent->PatchDataArraySelection
          ->ArrayExists(surfaceMeshName.c_str()))
        {
        this->Parent->PatchDataArraySelection
            ->DisableArray(surfaceMeshName.c_str());
        }
      this->SurfaceMeshSelectionStatus
          = this->Parent->GetPatchArrayStatus(surfaceMeshName.c_str());

      if (this->Parent->OutputProcessorPatches
          == vtkOFFReader::PROCESSOR_PATCHES_AGGREGATE)
        {
        const vtkStdString selectionName(this->RegionPrefix()
            + vtkOFFReaderPrivate::ProcessorPatchesIdentifier);
        if (!this->Parent->PatchDataArraySelection->ArrayExists(
            selectionName.c_str()))
          {
          this->Parent->PatchDataArraySelection->DisableArray(
              selectionName.c_str());
          }
        }

      // when in decommposed case (appended) mode, memorize processor
      // boundaries for quick lookup afterwards as processor boundary
      // faces has to be included in surface mesh appropriately
      const bool doProcBoundaries = this->SurfaceMeshSelectionStatus &&
          this->ProcessorName != "" && this->Parent->OutputProcessorPatches
          == vtkOFFReader::PROCESSOR_PATCHES_OFF;

      // iterate through each entry in the boundary file
      int allBoundariesNextStartFace = 0;
      this->BoundaryDict.resize(boundaryDict->size());
      for (size_t i = 0; i < boundaryDict->size(); i++)
        {
        vtkFoamEntry *boundaryEntryI = boundaryDict->operator[](i);
        const vtkFoamEntry *nFacesEntry = boundaryEntryI->Dictionary().Lookup("nFaces");
        if (nFacesEntry == NULL)
          {
          vtkErrorMacro(<< "nFaces entry not found in boundary entry "
              << boundaryEntryI->GetKeyword().c_str());
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          delete boundaryDict;
          return 0;
          }
        const int nFaces = nFacesEntry->ToInt();

        // extract name of the current patch for insertion
        const vtkStdString &boundaryNameI = boundaryEntryI->GetKeyword();

        // create BoundaryDict entry
        vtkFoamBoundaryEntry &BoundaryEntryI = this->BoundaryDict[i];
        BoundaryEntryI.NFaces = nFaces;
        BoundaryEntryI.BoundaryName = boundaryNameI;
        const vtkFoamEntry *startFaceEntry = boundaryEntryI->Dictionary().Lookup("startFace");
        if (startFaceEntry == NULL)
          {
          vtkErrorMacro(<< "startFace entry not found in boundary entry "
              << boundaryNameI.c_str());
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          delete boundaryDict;
          return 0;
          }
        BoundaryEntryI.StartFace = startFaceEntry->ToInt();
        const vtkFoamEntry *typeEntry = boundaryEntryI->Dictionary().Lookup("type");
        if (typeEntry == NULL)
          {
          vtkErrorMacro(<< "type entry not found in boundary entry "
              << boundaryNameI.c_str());
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          delete boundaryDict;
          return 0;
          }
        BoundaryEntryI.AllBoundariesStartFace = allBoundariesNextStartFace;
        const vtkStdString typeNameI(typeEntry->ToWord());
        // if the basic type of the patch is one of the followings the
        // point-filtered values at patches are overridden by patch values
        if (typeNameI == "patch" || typeNameI == "wall" || typeNameI == "ggi"
            || typeNameI == "cyclicGgi" || typeNameI == "overlapGgi")
          {
          BoundaryEntryI.BoundaryType = vtkFoamBoundaryEntry::PHYSICAL;
          allBoundariesNextStartFace += nFaces;
          }
        else if (typeNameI == "processor")
          {
          BoundaryEntryI.BoundaryType = vtkFoamBoundaryEntry::PROCESSOR;
          allBoundariesNextStartFace += nFaces;

          const vtkFoamEntry *procNoEntry
              = boundaryEntryI->Dictionary().Lookup("myProcNo");
          if (!procNoEntry)
            {
            vtkErrorMacro(<< "myProcNo entry not found in boundary entry "
                << boundaryNameI.c_str());
            this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
            delete boundaryDict;
            return 0;
            }
          BoundaryEntryI.MyProcNo = procNoEntry->ToInt();
          procNoEntry = boundaryEntryI->Dictionary().Lookup("neighbProcNo");
          if (!procNoEntry)
            {
            vtkErrorMacro(<< "neighbProcNo entry not found in boundary entry "
                << boundaryEntryI->GetKeyword().c_str());
            this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
            delete boundaryDict;
            return 0;
            }
          BoundaryEntryI.NeighbProcNo = procNoEntry->ToInt();

          if (BoundaryEntryI.MyProcNo > BoundaryEntryI.NeighbProcNo)
            {
            if (BoundaryEntryI.StartFace
                < this->BoundaryDict.LowestUpperProcFaceNo)
              {
              this->BoundaryDict.LowestUpperProcFaceNo
                  = BoundaryEntryI.StartFace;
              }
            const size_t rangeSize = this->BoundaryDict.UpperProcRanges.size();
            if (rangeSize > 1
                && this->BoundaryDict.UpperProcRanges[rangeSize - 1]
                == BoundaryEntryI.StartFace)
              {
              this->BoundaryDict.UpperProcRanges[rangeSize - 1]
                = BoundaryEntryI.StartFace + nFaces;
              }
            else
              {
              this->BoundaryDict.UpperProcRanges.push_back(
                  BoundaryEntryI.StartFace);
              this->BoundaryDict.UpperProcRanges.push_back(
                  BoundaryEntryI.StartFace + nFaces);
              }
            }

          if (doProcBoundaries)
            {
            this->BoundaryDict.ProcBoundaries.push_back(i);
            }
          }
        else
          {
          BoundaryEntryI.BoundaryType = vtkFoamBoundaryEntry::GEOMETRICAL;
          }
        BoundaryEntryI.IsActive = false;

        // always hide processor patches for decomposed cases (appended) to keep
        // vtkAppendCompositeDataLeaves happy
        if (this->Parent->OutputProcessorPatches
            == vtkOFFReader::PROCESSOR_PATCHES_OFF
            && BoundaryEntryI.BoundaryType == vtkFoamBoundaryEntry::PROCESSOR)
          {
          continue;
          }
        const vtkStdString selectionName(this->RegionPrefix()
            + (BoundaryEntryI.BoundaryType == vtkFoamBoundaryEntry::PROCESSOR
            && this->Parent->OutputProcessorPatches
            == vtkOFFReader::PROCESSOR_PATCHES_AGGREGATE
            ? vtkOFFReaderPrivate::ProcessorPatchesIdentifier
            : boundaryNameI));
        if (this->Parent->PatchDataArraySelection->
            ArrayExists(selectionName.c_str()))
          {
          // Mark boundary if selected for display
          if (this->Parent->GetPatchArrayStatus(selectionName.c_str()))
            {
            BoundaryEntryI.IsActive = true;
            }
          }
        else
          {
          // add patch to list with selection status turned off:
          // the patch is added to list even if its size is zero
          this->Parent->PatchDataArraySelection->DisableArray(selectionName.c_str());
          }
        }

      delete boundaryDict;
      }
    }

  // Add scalars and vectors to metadata
  vtkStdString timePath(this->CurrentTimePath());
  // do not do "RemoveAllArrays()" to accumulate array selections
  // this->CellDataArraySelection->RemoveAllArrays();
  this->VolFieldFiles->Initialize();
  this->SurfaceFieldFiles->Initialize();
  this->PointFieldFiles->Initialize();
  vtkStringArray *cellObjectNames = vtkStringArray::New();
  vtkStringArray *surfaceObjectNames = vtkStringArray::New();
  vtkStringArray *pointObjectNames = vtkStringArray::New();
  this->GetFieldNames(timePath + this->RegionPath(), false, cellObjectNames,
      surfaceObjectNames, pointObjectNames);

  this->LagrangianFieldFiles->Initialize();
  if (listNextTimeStep)
    {
    this->Parent->LagrangianPaths->Initialize();
    }
  vtkStringArray *lagrangianObjectNames = vtkStringArray::New();
  this->LocateLagrangianClouds(lagrangianObjectNames, timePath);

  // if the requested timestep is 0 then we also look at the next
  // timestep to add extra objects that don't exist at timestep 0 into
  // selection lists. Note the ObjectNames array will be recreated in
  // RequestData() so we don't have to worry about duplicated fields.
  if (listNextTimeStep && this->TimeValues->GetNumberOfTuples() >= 2
      && this->TimeStep == 0)
    {
    const vtkStdString timePath2(this->TimePath(1));
    this->GetFieldNames(timePath2 + this->RegionPath(), false, cellObjectNames,
        surfaceObjectNames, pointObjectNames);
    // if lagrangian clouds were not found at timestep 0
    if (this->Parent->LagrangianPaths->GetNumberOfTuples() == 0)
      {
      this->LocateLagrangianClouds(lagrangianObjectNames, timePath2);
      }
    }

  // sort array names
  this->SortFieldFiles(cellSelectionNames, this->VolFieldFiles, cellObjectNames);
  this->SortFieldFiles(surfaceSelectionNames, this->SurfaceFieldFiles,
      surfaceObjectNames);
  this->SortFieldFiles(pointSelectionNames, this->PointFieldFiles,
      pointObjectNames);
  this->SortFieldFiles(lagrangianSelectionNames, this->LagrangianFieldFiles,
      lagrangianObjectNames);

  return 1;
}

//-----------------------------------------------------------------------------
// list time directories according to controlDict
bool vtkOFFReaderPrivate::ListTimeDirectoriesByControlDict(
    vtkFoamDict* dictPtr)
{
  vtkFoamDict& dict = *dictPtr;

  const vtkFoamEntry *startTimeEntry = dict.Lookup("startTime");
  if (startTimeEntry == NULL)
    {
    vtkErrorMacro(<< "startTime entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  // using double to precisely handle time values
  const double startTime = startTimeEntry->ToDouble();

  const vtkFoamEntry *endTimeEntry = dict.Lookup("endTime");
  if (endTimeEntry == NULL)
    {
    vtkErrorMacro(<< "endTime entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  const double endTime = endTimeEntry->ToDouble();

  const vtkFoamEntry *deltaTEntry = dict.Lookup("deltaT");
  if (deltaTEntry == NULL)
    {
    vtkErrorMacro(<< "deltaT entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  const double deltaT = deltaTEntry->ToDouble();

  const vtkFoamEntry *writeIntervalEntry = dict.Lookup("writeInterval");
  if (writeIntervalEntry == NULL)
    {
    vtkErrorMacro(<< "writeInterval entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  const double writeInterval = writeIntervalEntry->ToDouble();

  const vtkFoamEntry *timeFormatEntry = dict.Lookup("timeFormat");
  if (timeFormatEntry == NULL)
    {
    vtkErrorMacro(<< "timeFormat entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  const vtkStdString timeFormat(timeFormatEntry->ToWord());

  const vtkFoamEntry *timePrecisionEntry = dict.Lookup("timePrecision");
  const int timePrecision // default is 6
      = (timePrecisionEntry != NULL ? timePrecisionEntry->ToInt() : 6);

  // calculate the time step increment based on type of run
  const vtkFoamEntry *writeControlEntry = dict.Lookup("writeControl");
  if (writeControlEntry == NULL)
    {
    vtkErrorMacro(<< "writeControl entry not found in controlDict");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }
  const vtkStdString writeControl(writeControlEntry->ToWord());
  double timeStepIncrement;
  if (writeControl == "timeStep")
    {
    timeStepIncrement = writeInterval * deltaT;
    }
  else if (writeControl == "runTime" || writeControl == "adjustableRunTime")
    {
    timeStepIncrement = writeInterval;
    }
  else
    {
    vtkErrorMacro(<<"Time step can't be determined because writeControl is"
        " set to " << writeControl.c_str());
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return false;
    }

  // calculate how many timesteps there should be
  const double tempResult = (endTime - startTime) / timeStepIncrement;
  // +0.5 to round up
  const int tempNumTimeSteps = static_cast<int>(tempResult + 0.5) + 1;

  // make sure time step dir exists
  vtkstd::vector<double> tempSteps;
  vtkDirectory *test = vtkDirectory::New();
  this->TimeValues->Initialize();
  this->TimeNames->Initialize();

  // determine time name based on Foam::Time::timeName()
  // cf. src/OpenFOAM/db/Time/Time.C
  vtksys_ios::ostringstream parser;
#ifdef _MSC_VER
  bool correctExponent = true;
#endif
  if (timeFormat == "general")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(0), vtksys_ios::ios_base::floatfield);
    }
  else if (timeFormat == "fixed")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(vtksys_ios::ios_base::fixed),
        vtksys_ios::ios_base::floatfield);
#ifdef _MSC_VER
    correctExponent = false;
#endif
    }
  else if (timeFormat == "scientific")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(vtksys_ios::ios_base::scientific),
        vtksys_ios::ios_base::floatfield);
    }
  else
    {
    vtkWarningMacro("Warning: unsupported time format. Assuming general.");
    parser.setf(vtksys_ios::ios_base::fmtflags(0), vtksys_ios::ios_base::floatfield);
    }
  parser.precision(timePrecision);

  for (int i = 0; i < tempNumTimeSteps; i++)
    {
    parser.str("");
    const double tempStep = i * timeStepIncrement + startTime;
    parser << tempStep; // stringstream doesn't require ends
#ifdef _MSC_VER
    // workaround for format difference in MSVC++:
    // remove an extra 0 from exponent
    if(correctExponent)
      {
      vtkStdString tempStr(parser.str());
      vtkStdString::size_type pos = tempStr.find('e');
      if(pos != vtkStdString::npos && tempStr.length() >= pos + 3
          && tempStr[pos + 2] == '0')
        {
        tempStr.erase(pos + 2, 1);
        parser.str(tempStr);
        }
      }
#endif
    // Add the time steps that actually exist to steps
    // allows the run to be stopped short of controlDict spec
    // allows for removal of timesteps
    if (test->Open((this->CasePath + parser.str()).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(parser.str());
      }
    // necessary for reading the case/0 directory whatever the timeFormat is
    // based on Foam::Time::operator++() cf. src/OpenFOAM/db/Time/Time.C
    else if ((fabs(tempStep) < 1.0e-14L) // 10*SMALL
        && test->Open((this->CasePath + vtkStdString("0")).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(vtkStdString("0"));
      }
    }
  test->Delete();
  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  if (this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if (test->Open((this->CasePath + "constant").c_str()))
      {
      parser.str("");
      parser << startTime;
      this->TimeValues->InsertNextValue(startTime);
      this->TimeValues->Squeeze();
      this->TimeNames->InsertNextValue(parser.str());
      this->TimeNames->Squeeze();
      }
    test->Delete();
    }
  return true;
}

//-----------------------------------------------------------------------------
// list time directories by searching all valid time instances in a
// case directory
bool vtkOFFReaderPrivate::ListTimeDirectoriesByInstances()
{
  // open the case directory
  vtkDirectory* test = vtkDirectory::New();
  if (!test->Open(this->CasePath.c_str()))
    {
    test->Delete();
    vtkErrorMacro(<< "Can't open directory " << this->CasePath.c_str());
    this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    return false;
    }

  // search all the directories in the case directory and detect
  // directories with names convertible to numbers
  this->TimeValues->Initialize();
  this->TimeNames->Initialize();
  const int nFiles = test->GetNumberOfFiles();
  for (int i = 0; i < nFiles; i++)
    {
    const vtkStdString dir = test->GetFile(i);
    if (test->FileIsDirectory(dir.c_str()))
      {
      // check if the name is convertible to a number
      bool isTimeDir = true;
      for (size_t j = 0; j < dir.length(); j++)
        {
        const char c = dir[j];
        if (!isdigit(c) && c != '+' && c != '-' && c != '.' && c != 'e' && c
            != 'E')
          {
          isTimeDir = false;
          break;
          }
        }
      if (!isTimeDir)
        {
        continue;
        }

      // convert to a number
#if VTK_FOAMFILE_LOCALE_WORKAROUND
      vtksys_ios::istringstream timeStream(dir);
      double timeValue;
      timeStream >> timeValue;

      // check if the value really was converted to a number
      if (timeStream.fail() || !timeStream.eof())
        {
        continue;
        }
#else
      char *endptr;
      double timeValue = strtod(dir.c_str(), &endptr);
      // check if the value really was converted to a number
      if (timeValue == 0.0 && endptr == dir.c_str())
        {
        continue;
        }
#endif

      // add to the instance list
      this->TimeValues->InsertNextValue(timeValue);
      this->TimeNames->InsertNextValue(dir);
      }
    }
  test->Delete();

  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  if (this->TimeValues->GetNumberOfTuples() > 1)
    {
    // sort the detected time directories
    vtkSortDataArray::Sort(this->TimeValues, this->TimeNames);

    // if there are duplicated timeValues found, remove duplicates
    // (e.g. "0" and "0.000")
    for (int timeI = 1; timeI < this->TimeValues->GetNumberOfTuples();)
      {
      // compare by exact match
      if (this->TimeValues->GetValue(timeI - 1)
          == this->TimeValues->GetValue(timeI))
        {
        vtkWarningMacro(<<"Different time directories with the same time value "
            << this->TimeNames->GetValue(timeI - 1).c_str() << " and "
            << this->TimeNames->GetValue(timeI).c_str() << " found. "
            << this->TimeNames->GetValue(timeI).c_str() << " will be ignored.");
        this->TimeValues->RemoveTuple(timeI);
        // vtkStringArray does not have RemoveTuple()
        for (int timeJ = timeI + 1; timeJ
            < this->TimeNames->GetNumberOfTuples(); timeJ++)
          {
          this->TimeNames->SetValue(timeJ - 1, this->TimeNames->GetValue(timeJ));
          }
        this->TimeNames->Resize(this->TimeNames->GetNumberOfTuples() - 1);
        }
      else
        {
        timeI++;
        }
      }
    }
  else if (this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if (test->Open((this->CasePath + "constant").c_str()))
      {
      this->TimeValues->InsertNextValue(0.0);
      this->TimeValues->Squeeze();
      this->TimeNames->InsertNextValue("0");
      this->TimeNames->Squeeze();
      }
    test->Delete();
    }

  return true;
}

//-----------------------------------------------------------------------------
// gather the necessary information to create a path to the data
bool vtkOFFReaderPrivate::MakeInformationVector(
    const vtkStdString &casePath, const vtkStdString &controlDictPath,
    const vtkStdString &procName, vtkOFFReader *parent)
{
  this->CasePath = casePath;
  this->ProcessorName = procName;
  this->Parent = parent;

  // list timesteps (skip parsing controlDict entirely if
  // ListTimeStepsByControlDict is set to 0)
  bool ret = false; // tentatively set to false to suppress warning by older compilers
  if (this->Parent->GetListTimeStepsByControlDict())
    {
    vtkFoamIOobject io(this->CasePath,
        this->Parent->GetIsSinglePrecisionBinary() != 0);

    // open and check if controlDict is readable
    if (!io.Open(controlDictPath))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
      return false;
      }
    vtkFoamDict dict;
    if (!dict.Read(io))
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }
    if (dict.GetType() != vtkFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<<"The file type of " << io.GetFileName().c_str()
          << " is not a dictionary");
      this->Parent->SetErrorCode(vtkErrorCode::UnrecognizedFileTypeError);
      return false;
      }

    const vtkFoamEntry *writeControlEntry = dict.Lookup("writeControl");
    if (writeControlEntry == NULL)
      {
      vtkErrorMacro(<< "writeControl entry not found in "
          << io.GetFileName().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }
    const vtkStdString writeControl(writeControlEntry->ToWord());

    // empty if not found
    const vtkFoamEntry *adjustTimeStepEntry = dict.Lookup("adjustTimeStep");
    const vtkStdString
        adjustTimeStep = adjustTimeStepEntry == NULL ? vtkStdString()
            : adjustTimeStepEntry->ToWord();

    // list time directories according to controlDict if (adjustTimeStep
    // writeControl) == (off, timeStep) or (on, adjustableRunTime); list
    // by time instances in the case directory otherwise (different behavior
    // from paraFoam)
    // valid switching words cf. src/OpenFOAM/db/Switch/Switch.C
    if ((((adjustTimeStep == "off" || adjustTimeStep == "no" || adjustTimeStep
        == "n" || adjustTimeStep == "false" || adjustTimeStep == "")
        && writeControl == "timeStep") || ((adjustTimeStep == "on"
        || adjustTimeStep == "yes" || adjustTimeStep == "y" || adjustTimeStep
        == "true") && writeControl == "adjustableRunTime")))
      {
      ret = this->ListTimeDirectoriesByControlDict(&dict);
      }
    else
      {
      ret = this->ListTimeDirectoriesByInstances();
      }
    }
  else
    {
    ret = this->ListTimeDirectoriesByInstances();
    }

  if (!ret)
    {
    return ret;
    }

  // does not seem to be required even if number of timesteps reduced
  // upon refresh since ParaView rewinds TimeStep to 0, but for precaution
  if (this->TimeValues->GetNumberOfTuples() > 0)
    {
    if (this->TimeStep >= this->TimeValues->GetNumberOfTuples())
      {
      this->SetTimeStep(this->TimeValues->GetNumberOfTuples() - 1);
      }
    }
  else
    {
    this->SetTimeStep(0);
    }

  this->PopulatePolyMeshDirArrays();
  return ret;
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::AppendMeshDirToArray(
    vtkStringArray* polyMeshDir, const vtkStdString &path,
    const vtkStdString &fileName, const int timeI)
{
  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkStdString filePath(path + fileName);
  if (io.Open(filePath) || io.Open(filePath + ".gz"))
    {
    io.Close();
    // set points/faces location to current timesteps value
    polyMeshDir->SetValue(timeI, this->TimeNames->GetValue(timeI));
    }
  else
    {
    if (timeI != 0)
      {
      // set points/faces location to previous timesteps value
      polyMeshDir->SetValue(timeI, polyMeshDir->GetValue(timeI - 1));
      }
    else
      {
      filePath = this->CasePath + "constant" + this->RegionPath()
          + "/polyMesh/" + fileName;
      if(io.Open(filePath) || io.Open(filePath + ".gz"))
        {
        // set points/faces to constant
        polyMeshDir->SetValue(timeI, "constant");
        }
      }
    }
}

//-----------------------------------------------------------------------------
// create a Lookup Table containing the location of the points
// and faces files for each time steps mesh
void vtkOFFReaderPrivate::PopulatePolyMeshDirArrays()
{
  // intialize size to number of timesteps
  const int nSteps = this->TimeValues->GetNumberOfTuples();
  this->PolyMeshPointsDir->SetNumberOfValues(nSteps);
  this->PolyMeshFacesDir->SetNumberOfValues(nSteps);

  // loop through each timestep
  for (int i = 0; i < nSteps; i++)
    {
    // create the path to the timestep
    vtkStdString polyMeshPath = this->TimeRegionPath(i) + "/polyMesh/";
    AppendMeshDirToArray(this->PolyMeshPointsDir, polyMeshPath, "points", i);
    AppendMeshDirToArray(this->PolyMeshFacesDir, polyMeshPath, "faces", i);
    }
  return;
}

//-----------------------------------------------------------------------------
// read the points file into a vtkFloatArray
vtkFloatArray* vtkOFFReaderPrivate::ReadPointsFile()
{
  if (this->PolyMeshPointsDir->GetValue(this->TimeStep) == "")
    {
    vtkErrorMacro(<<"Cannot find path to points file for region \""
        << this->RegionName.c_str() << "\" in case " << this->CasePath.c_str());
    this->Parent->SetErrorCode(vtkErrorCode::UnknownError);
    return NULL;
    }

  // path to points file
  const vtkStdString pointPath =
      this->CurrentTimeRegionMeshPath(this->PolyMeshPointsDir) + "points";

  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(pointPath) || io.Open(pointPath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str());
    this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    return NULL;
    }

  vtkFoamEntryValue dict(NULL);
  try
    {
    dict.ReadNonuniformList<vtkFoamToken::VECTORLIST,
    vtkFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, false> >(io);
    }
  catch(vtkFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << e.c_str());
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return NULL;
    }

  vtkFloatArray *pointArray = static_cast<vtkFloatArray *>(dict.Ptr());

  // set the number of points
  this->NumPoints = pointArray->GetNumberOfTuples();

  return pointArray;
}

//-----------------------------------------------------------------------------
// read the faces into a vtkFoamIntVectorVector
vtkFoamIntVectorVector * vtkOFFReaderPrivate::ReadFacesFile(
    const vtkStdString &facePathIn)
{
  const vtkStdString facePath(facePathIn + "faces");

  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(facePath) || io.Open(facePath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str() << ". If you are trying to read a parallel "
        "decomposed case, set Case Type to Decomposed Case.");
    this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    return NULL;
    }

  vtkFoamEntryValue dict(NULL);
  try
    {
    if (io.GetClassName() == "faceCompactList")
      {
      dict.ReadCompactIOLabelList(io);
      }
    else
      {
      dict.ReadLabelListList(io);
      }
    }
  catch(vtkFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << e.c_str());
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return NULL;
    }
  return static_cast<vtkFoamIntVectorVector *>(dict.Ptr());
}

//-----------------------------------------------------------------------------
// read the owner and neighbor file and create cellFaces
vtkFoamIntVectorVector * vtkOFFReaderPrivate::ReadOwnerNeighborFiles(
    const vtkStdString &ownerNeighborPath, vtkFoamIntVectorVector *facePoints)
{
  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkStdString ownerPath(ownerNeighborPath + "owner");
  if (io.Open(ownerPath) || io.Open(ownerPath + ".gz"))
    {
    vtkFoamEntryValue ownerDict(NULL);
    try
      {
      ownerDict.ReadNonuniformList<vtkFoamToken::LABELLIST,
      vtkFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }

    io.Close();

    const vtkStdString neighborPath(ownerNeighborPath + "neighbour");
    if (!(io.Open(neighborPath) || io.Open(neighborPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
      return NULL;
      }

    vtkFoamEntryValue neighborDict(NULL);
    try
      {
      neighborDict.ReadNonuniformList<vtkFoamToken::LABELLIST,
      vtkFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }

    this->FaceOwner = static_cast<vtkIntArray *>(ownerDict.Ptr());
    vtkIntArray &faceOwner = *this->FaceOwner;
    vtkIntArray &faceNeighbor = neighborDict.LabelList();

    const int nFaces = faceOwner.GetNumberOfTuples();
    const int nNeiFaces = faceNeighbor.GetNumberOfTuples();

    if (nFaces < nNeiFaces)
      {
      vtkErrorMacro(<<"Number of owner faces " << nFaces
          << " must be equal or larger than number of neighbor faces "
          << nNeiFaces);
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }

    // the number of owner faces can be smaller than the number of
    // faces in sliding mesh cases
    if (nFaces > facePoints->GetNumberOfElements())
      {
      vtkWarningMacro(<<"Number of faces in faces "
          << facePoints->GetNumberOfElements()
          << " must be equal or larger than number of owner faces "
          << nFaces);
      return NULL;
      }

    // add the face numbers to the correct cell cf. Terry's code and
    // src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCells.C
    // find the number of cells
    int nCells = -1;
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      // we do need to take neighbor faces into account since all the
      // surrounding faces of a cell can be neighbors for a valid mesh
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if (nCells < neighborCell) // max(nCells, faceNeighbor[i])
        {
        nCells = neighborCell;
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      }
    nCells++;

    if (nCells == 0)
      {
      vtkWarningMacro(<<"The mesh contains no cells");
      }

    // set the number of cells
    this->NumCells = nCells;

    // create cellFaces with the length of the body undetermined
    vtkFoamIntVectorVector *cells = new vtkFoamIntVectorVector(nCells, 1);

    // count number of faces for each cell
    int *cfiPtr = cells->GetIndices()->GetPointer(0);
    for (int cellI = 0; cellI <= nCells; cellI++)
      {
      cfiPtr[cellI] = 0;
      }
    int nTotalCellFaces = 0;
    cfiPtr++; // offset +1
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if (neighborCell >= 0)
        {
        cfiPtr[neighborCell]++;
        nTotalCellFaces++;
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      }
    cfiPtr--; // revert offset +1

    // allocate cellFaces. To reduce the numbers of new/delete operations we
    // allocate memory space for all faces linearly
    cells->ResizeBody(nTotalCellFaces);

    // accumulate the number of cellFaces to create cellFaces indices
    // and copy them to a temporary array
    vtkIntArray *tmpFaceIndices = vtkIntArray::New();
    tmpFaceIndices->SetNumberOfValues(nCells + 1);
    int *tfiPtr = tmpFaceIndices->GetPointer(0);
    tfiPtr[0] = 0;
    for (int cellI = 1; cellI <= nCells; cellI++)
      {
      tfiPtr[cellI] = (cfiPtr[cellI] += cfiPtr[cellI - 1]);
      }

    // add face numbers to cell-faces list
    vtkIntArray *cellFacesList = cells->GetBody();
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[ownerCell]++, faceI);
        }
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if (neighborCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[neighborCell]++, faceI);
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[ownerCell]++, faceI);
        }
      }
    tmpFaceIndices->Delete();

    return cells;
    }
  else // if owner does not exist look for cells
    {
    vtkStdString cellsPath(ownerNeighborPath + "cells");
    if (!(io.Open(cellsPath) || io.Open(cellsPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
      return NULL;
      }
    vtkFoamEntryValue cellsDict(NULL);
    try
      {
      cellsDict.ReadLabelListList(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }

    vtkFoamIntVectorVector *cells =
        static_cast<vtkFoamIntVectorVector *>(cellsDict.Ptr());
    this->NumCells = cells->GetNumberOfElements();
    const int nFaces = facePoints->GetNumberOfElements();

    // create face owner list
    this->FaceOwner = vtkIntArray::New();
    this->FaceOwner->SetNumberOfTuples(nFaces);
    for (int faceI = 0; faceI < nFaces; faceI++)
      {
      this->FaceOwner->SetValue(faceI, -1);
      }
    for (int cellI = 0; cellI < this->NumCells; cellI++)
      {
      const int nCellFaces = cells->GetSize(cellI);
      const int *cellFaces = cells->operator[](cellI);
      for (int faceI = 0; faceI < nCellFaces; faceI++)
        {
        const int f = cellFaces[faceI];
        if (f < 0 || f >= nFaces) // make sure the face number is valid
          {
          vtkErrorMacro("Face number " << f << " in cell " << cellI
              << " exceeds the number of faces " << nFaces);
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          this->FaceOwner->Delete();
          this->FaceOwner = NULL;
          delete cells;
          return NULL;
          }
        const int owner = this->FaceOwner->GetValue(f);
        if (owner == -1 || owner > cellI)
          {
          this->FaceOwner->SetValue(f, cellI);
          }
        }
      }
    // check for unused faces
    for (int faceI = 0; faceI < nFaces; faceI++)
      {
      if (this->FaceOwner->GetValue(faceI) == -1)
        {
        vtkErrorMacro(<<"Face " << faceI << " is not used");
        this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
        this->FaceOwner->Delete();
        this->FaceOwner = NULL;
        delete cells;
        return NULL;
        }
      }

    return cells;
    }
}

//-----------------------------------------------------------------------------
bool vtkOFFReaderPrivate::CheckFacePoints(
    vtkFoamIntVectorVector *facePoints)
{
  const int nFaces = facePoints->GetNumberOfElements();

  for (int faceI = 0; faceI < nFaces; faceI++)
    {
    const int nPoints = facePoints->GetSize(faceI);
    const int *pointList = facePoints->operator[](faceI);
    if (nPoints < 3)
      {
      vtkErrorMacro(<< "Face " << faceI << " has only " << nPoints
          << " points which is not enough to constitute a face"
          " (a face must have at least 3 points)");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int p = pointList[pointI];
      if (p < 0 || p >= this->NumPoints)
        {
        vtkErrorMacro(<< "The point number " << p << " at face number " << faceI
            << " is out of range for " << this->NumPoints << " points");
        this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
        return false;
        }
      }
    }
  return true;
}

//-----------------------------------------------------------------------------
// determine cell shape and insert the cell into the mesh
// hexahedron, prism, pyramid, tetrahedron and decompose polyhedron
void vtkOFFReaderPrivate::InsertCellsToGrid(
    vtkUnstructuredGrid* internalMesh,
    const vtkFoamIntVectorVector *cellsFaces,
    const vtkFoamIntVectorVector *facesPoints, vtkFloatArray *pointArray,
    vtkIdTypeArray *additionalCells, vtkIntArray *cellList)
{
  const int maxNPoints = 256; // assume max number of points per cell
  vtkIdList* cellPoints = vtkIdList::New();
  cellPoints->SetNumberOfIds(maxNPoints);
#if VTK_FOAMFILE_USE_VTK_POLYHEDRON
  // assume max number of nPoints per face + points per cell
  const int maxNPolyPoints = 1024;
  vtkIdList* polyPoints = vtkIdList::New();
  polyPoints->SetNumberOfIds(maxNPolyPoints);
#endif

  const int nCells = (cellList == NULL ? this->NumCells
      : cellList->GetNumberOfTuples());
  int nAdditionalPoints = 0;
  this->NumTotalAdditionalCells = 0;

  // alias
  const vtkFoamIntVectorVector& facePoints = *facesPoints;

  for (int cellI = 0; cellI < nCells; cellI++)
    {
    int cellId;
    if (cellList == NULL)
      {
      cellId = cellI;
      }
    else
      {
      cellId = cellList->GetValue(cellI);
      if (cellId >= this->NumCells)
        {
        vtkWarningMacro(<<"cellLabels id " << cellId
            << " exceeds the number of cells " << nCells
            << ". Inserting an empty cell.");
        internalMesh->InsertNextCell(VTK_EMPTY_CELL, 0,
            cellPoints->GetPointer(0));
        continue;
        }
      }
    const int *cellFaces = cellsFaces->operator[](cellId);
    const int nCellFaces = cellsFaces->GetSize(cellId);

    // determine type of the cell
    // cf. src/OpenFOAM/meshes/meshShapes/cellMatcher/{hex|prism|pyr|tet}-
    // Matcher.C
    int cellType = VTK_CONVEX_POINT_SET;
    if (nCellFaces == 6)
      {
      int j = 0;
      for (; j < nCellFaces; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) != 4)
          {
          break;
          }
        }
      if (j == nCellFaces)
        {
        cellType = VTK_HEXAHEDRON;
        }
      }
    else if (nCellFaces == 5)
      {
      int nTris = 0, nQuads = 0;
      for (int j = 0; j < nCellFaces; j++)
        {
        const int nPoints = facePoints.GetSize(cellFaces[j]);
        if (nPoints == 3)
          {
          nTris++;
          }
        else if (nPoints == 4)
          {
          nQuads++;
          }
        else
          {
          break;
          }
        }
      if (nTris == 2 && nQuads == 3)
        {
        cellType = VTK_WEDGE;
        }
      else if (nTris == 4 && nQuads == 1)
        {
        cellType = VTK_PYRAMID;
        }
      }
    else if (nCellFaces == 4)
      {
      int j = 0;
      for (; j < nCellFaces; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) != 3)
          {
          break;
          }
        }
      if (j == nCellFaces)
        {
        cellType = VTK_TETRA;
        }
      }

    // not a Hex/Wedge/Pyramid/Tetra
    if (cellType == VTK_CONVEX_POINT_SET)
      {
      int nPoints = 0;
      for (int j = 0; j < nCellFaces; j++)
        {
        nPoints += facePoints.GetSize(cellFaces[j]);
        }
      if (nPoints == 0)
        {
        cellType = VTK_EMPTY_CELL;
        }
      }

    // Cell shape constructor based on the one implementd by Terry
    // Jordan, with lots of improvements. Not as elegant as the one in
    // OpenFOAM but it's simple and works reasonably fast.

    // OFhex | vtkHexahedron
    if (cellType == VTK_HEXAHEDRON)
      {
      // get first face in correct order
      const int cellBaseFaceId = cellFaces[0];
      const int *face0Points = facePoints[cellBaseFaceId];

      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for (int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[3 - j]);
          }
        }
      else
        {
        // add base face to cell points
        for (int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      int dupPoint = -1;
      for (int faceI = 1; faceI < 5; faceI++) // skip face 0 and 5
        {
        const int cellFaceI = cellFaces[faceI];
        const int *faceIPoints = facePoints[cellFaceI];
        int foundDup = -1, pointI = 0;
        for (; pointI < 4; pointI++) // each point
          {
          const int faceIPointI = faceIPoints[pointI];
          // matching two points in base face is enough to find a
          // duplicated point since neighboring faces share two
          // neighboring points (i. e. an edge)
          if (baseFacePoint0 == faceIPointI)
            {
            foundDup = 0;
            break;
            }
          else if (baseFacePoint2 == faceIPointI)
            {
            foundDup = 2;
            break;
            }
          }
        if (foundDup >= 0)
          {
          // find the pivot point if still haven't
          if (pivotPoint == -1)
            {
            dupPoint = foundDup;

            const int faceINextPoint = faceIPoints[(pointI + 1) % 4];

            // if the next point of the faceI-th face matches the
            // previous point of the base face use the previous point
            // of the faceI-th face as the pivot point; or use the
            // next point otherwise
            if (faceINextPoint == (this->FaceOwner->GetValue(cellFaceI)
                == cellId ? cellPoints->GetId(1 + foundDup)
                : cellPoints->GetId(3 - foundDup)))
              {
              pivotPoint = faceIPoints[(3 + pointI) % 4];
              }
            else
              {
              pivotPoint = faceINextPoint;
              }

            if (cellOppositeFaceI >= 0)
              {
              break;
              }
            }
          }
        else
          {
          // if no duplicated point found, faceI is the opposite face
          cellOppositeFaceI = cellFaceI;

          if (pivotPoint >= 0)
            {
            break;
            }
          }
        }

      // if the opposite face is not found until face 4, face 5 is
      // always the opposite face
      if (cellOppositeFaceI == -1)
        {
        cellOppositeFaceI = cellFaces[5];
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for (; pivotPointI < 4; pivotPointI++)
        {
        if (oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      // shift the pivot point if the point corresponds to point 2
      // of the base face
      if (dupPoint == 2)
        {
        pivotPointI = (pivotPointI + 2) % 4;
        }
      // copy the face-point list of the opposite face to cell-point list
      int basePointI = 4;
      if (this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        for (int pointI = pivotPointI; pointI < 4; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 0; pointI < pivotPointI; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        for (int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 3; pointI > pivotPointI; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }

      // create the hex cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, 8, cellPoints->GetPointer(0));
      }

    // the cell construction is about the same as that of a hex, but
    // the point ordering have to be reversed!!
    else if (cellType == VTK_WEDGE)
      {
      // find the base face number
      int baseFaceId = 0;
      for (int j = 0; j < 5; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) == 3)
          {
          baseFaceId = j;
          break;
          }
        }

      // get first face in correct order
      const int cellBaseFaceId = cellFaces[baseFaceId];
      const int *face0Points = facePoints[cellBaseFaceId];

      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        for (int j = 0; j < 3; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      else
        {
        // if it is a neighbor face flip the points
        for (int j = 0; j < 3; j++)
          {
          // add base face to cell points
          cellPoints->SetId(j, face0Points[2 - j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      bool dupPoint2 = false;
      for (int faceI = 0; faceI < 5; faceI++)
        {
        if (faceI == baseFaceId)
          {
          continue;
          }
        const int cellFaceI = cellFaces[faceI];
        if (facePoints.GetSize(cellFaceI) == 3)
          {
          cellOppositeFaceI = cellFaceI;
          }
        // find the pivot point if still haven't
        else if (pivotPoint == -1)
          {
          const int *faceIPoints = facePoints[cellFaceI];
          bool found0Dup = false, found2Dup = false;
          int pointI = 0;
          for (; pointI < 4; pointI++) // each point
            {
            const int faceIPointI = faceIPoints[pointI];
            // matching two points in base face is enough to find a
            // duplicated point since neighboring faces share two
            // neighboring points (i. e. an edge)
            if (baseFacePoint0 == faceIPointI)
              {
              found0Dup = true;
              break;
              }
            else if (baseFacePoint2 == faceIPointI)
              {
              found2Dup = true;
              break;
              }
            }
          // the matching point must always be found so omit the check
          int baseFacePrevPoint, baseFaceNextPoint;
          if (found0Dup)
            {
            baseFacePrevPoint = cellPoints->GetId(2);
            baseFaceNextPoint = cellPoints->GetId(1);
            }
          else
            {
            baseFacePrevPoint = cellPoints->GetId(1);
            baseFaceNextPoint = cellPoints->GetId(0);
            dupPoint2 = true;
            }

          const int faceINextPoint = faceIPoints[(pointI + 1) % 4];
          const int faceIPrevPoint = faceIPoints[(3 + pointI) % 4];

          // if the next point of the faceI-th face matches the
          // previous point of the base face use the previous point of
          // the faceI-th face as the pivot point; or use the next
          // point otherwise
          if (faceINextPoint
              == (this->FaceOwner->GetValue(cellFaceI) == cellId ? baseFacePrevPoint
                  : baseFaceNextPoint))
            {
            pivotPoint = faceIPrevPoint;
            }
          else
            {
            pivotPoint = faceINextPoint;
            }
          }

        // break when both of opposite face and pivot point are found
        if (cellOppositeFaceI >= 0 && pivotPoint >= 0)
          {
          break;
          }
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for (; pivotPointI < 3; pivotPointI++)
        {
        if (oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      if (this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        if (dupPoint2)
          {
          pivotPointI = (pivotPointI + 2) % 3;
          }
        int basePointI = 3;
        for (int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 2; pointI > pivotPointI; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        // shift the pivot point if the point corresponds to point 2
        // of the base face
        if (dupPoint2)
          {
          pivotPointI = (1 + pivotPointI) % 3;
          }
        // copy the face-point list of the opposite face to cell-point list
        int basePointI = 3;
        for (int pointI = pivotPointI; pointI < 3; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 0; pointI < pivotPointI; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }

      // create the wedge cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, 6, cellPoints->GetPointer(0));
      }

    // OFpyramid | vtkPyramid || OFtet | vtkTetrahedron
    else if (cellType == VTK_PYRAMID || cellType == VTK_TETRA)
      {
      int baseFaceId = -1, nPoints;
      if (cellType == VTK_PYRAMID)
        {
        for (int j = 0; j < nCellFaces; j++)
          {
          if (facePoints.GetSize(cellFaces[j]) == 4)
            {
            baseFaceId = j;
            break;
            }
          }
        nPoints = 5;
        }
      else // VTK_TETRA
        {
        baseFaceId = 0;
        nPoints = 4;
        }

      // add first face to cell points
      const int cellBaseFaceId = cellFaces[baseFaceId];
      const int *baseFacePoints = facePoints[cellBaseFaceId];
      const vtkIdType nBaseFacePoints = facePoints.GetSize(cellBaseFaceId);
      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for (vtkIdType j = 0; j < nBaseFacePoints; j++)
          {
          cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
          }
        }
      else
        {
        for (vtkIdType j = 0; j < nBaseFacePoints; j++)
          {
          cellPoints->SetId(j, baseFacePoints[j]);
          }
        }

      // compare an adjacent face (any non base face is ok) point 1 to
      // base face points
      const int adjacentFaceId = (baseFaceId == 0) ? 1 : baseFaceId - 1;
      const int cellAdjacentFaceId = cellFaces[adjacentFaceId];
      const int *adjacentFacePoints = facePoints[cellAdjacentFaceId];
      const int adjacentFacePoint1 = adjacentFacePoints[1];
      bool foundDup = false;
      for (vtkIdType j = 0; j < nBaseFacePoints; j++)
        {
        // if point 1 of the adjacent face matches point j of the base face...
        if (cellPoints->GetId(j) == adjacentFacePoint1)
          {
          // if point 2 of the adjacent face matches the previous point
          // of the base face use point 0 of the adjacent face as the
          // pivot point; use point 2 otherwise
          cellPoints->SetId(
              nBaseFacePoints,
              (adjacentFacePoints[2]
                  == cellPoints->GetId((this->FaceOwner->GetValue(cellAdjacentFaceId)
                      == cellId ? (j + 1) : (nBaseFacePoints + j - 1))
                      % nBaseFacePoints)) ? adjacentFacePoints[0]
                  : adjacentFacePoints[2]);
          foundDup = true;
          break;
          }
        }
      // if point 1 of the adjacent face does not match any points of
      // the base face, it's the pivot point
      if (!foundDup)
        {
        cellPoints->SetId(nBaseFacePoints, adjacentFacePoint1);
        }

      // create the tetra cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, nPoints, cellPoints->GetPointer(0));
      }

    // erronous cells
    else if (cellType == VTK_EMPTY_CELL)
      {
      vtkWarningMacro("Warning: No points in cellId " << cellId);
      internalMesh->InsertNextCell(VTK_EMPTY_CELL, 0, cellPoints->GetPointer(0));
      }

    // OFpolyhedron || vtkConvexPointSet
    else
      {
      if (additionalCells != NULL) // decompose into tets and pyramids
        {
        // calculate cell centroid and insert it to point list
        this->AdditionalCellPoints->push_back(vtkIdList::New());
        vtkIdList *polyCellPoints = this->AdditionalCellPoints->back();
        float centroid[3];
        centroid[0] = centroid[1] = centroid[2] = 0.0F;
        for (int j = 0; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.GetSize(cellFacesJ);
          for (size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for (vtkIdType l = 0; l < static_cast<vtkIdType>(polyCellPoints->GetNumberOfIds()); l++)
              {
              if (polyCellPoints->GetId(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if (!foundDup)
              {
              polyCellPoints->InsertNextId(faceJPointK);
              const float *pointK = pointArray->GetPointer(3 * faceJPointK);
              centroid[0] += pointK[0];
              centroid[1] += pointK[1];
              centroid[2] += pointK[2];
              }
            }
          }
        polyCellPoints->Squeeze();
        const float weight = 1.0F
            / static_cast<float>(polyCellPoints->GetNumberOfIds());
        centroid[0] *= weight;
        centroid[1] *= weight;
        centroid[2] *= weight;
        pointArray->InsertNextTuple(centroid);

        // polyhedron decomposition.
        // a tweaked algorithm based on applications/utilities/postProcessing/
        // graphics/PVFoamReader/vtkFoam/vtkFoamAddInternalMesh.C
        bool insertDecomposedCell = true;
        int nAdditionalCells = 0;
        for (int j = 0; j < nCellFaces; j++)
          {
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const int nFaceJPoints = facePoints.GetSize(cellFacesJ);
          const int nTris = nFaceJPoints % 2;

          // check if the face has to be reversed
          bool reverseFace = false;
          if (cellFacesJ >= this->BoundaryDict.LowestUpperProcFaceNo)
            {
            vtkstd::vector<int> &ranges = this->BoundaryDict.UpperProcRanges;
            const size_t rangeSize = ranges.size();
            for (size_t rangeI = 0; rangeI < rangeSize; rangeI += 2)
              {
              if (ranges[rangeI] <= cellFacesJ
                  && cellFacesJ < ranges[rangeI + 1])
                {
                reverseFace = true;
                break;
                }
              }
            }

          int vertI = 2;

          // shift the start and end of the vertex loop if the
          // triangle of a decomposed face is going to be flat. Far
          // from perfect but better than nothing to avoid flat cells
          // which stops time integration of Stream Tracer especially
          // for split-hex unstructured meshes created by
          // e. g. autoRefineMesh
          if (nFaceJPoints >= 5 && nTris)
            {
            float *point0, *point1, *point2;
            point1 = pointArray->GetPointer(3 * faceJPoints[0]);
            if (reverseFace)
              {
              point0 = pointArray->GetPointer(3 * faceJPoints[1]);
              point2 = pointArray->GetPointer(3 * faceJPoints[2]);
              }
            else
              {
              point0
                  = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 1]);
              point2
                  = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 2]);
              }
            float vsizeSqr1 = 0.0F, vsizeSqr2 = 0.0F, dotProduct = 0.0F;
            for (int i = 0; i < 3; i++)
              {
              const float v1 = point1[i] - point0[i], v2 = point2[i]
                  - point0[i];
              vsizeSqr1 += v1 * v1;
              vsizeSqr2 += v2 * v2;
              dotProduct += v1 * v2;
              }
            // compare in squared representation to avoid using sqrt()
            if (dotProduct * (float) fabs(dotProduct) / (vsizeSqr1 * vsizeSqr2)
                < -1.0F + 1.0e-3F)
              {
              vertI = (reverseFace ? 3 : 1);
              }
            }

          if (reverseFace)
            {
            cellPoints->SetId(0, faceJPoints[vertI == 2 ? 0 : 1]);
            // if the number of vertices is odd there's a triangle
            if (nTris)
              {
              cellPoints->SetId(1, faceJPoints[vertI]);
              cellPoints->SetId(2, faceJPoints[vertI - 1]);
              cellPoints->SetId(3, this->NumPoints + nAdditionalPoints);

              if (insertDecomposedCell)
                {
                internalMesh->InsertNextCell(VTK_TETRA, 4,
                    cellPoints->GetPointer(0));
                insertDecomposedCell = false;
                }
              else
                {
                // set the 5th vertex number to -1 to distinguish a tetra cell
                cellPoints->SetId(4, -1);
                nAdditionalCells++;
                additionalCells->InsertNextTupleValue(
                    cellPoints->GetPointer(0));
                }

              ++vertI;
              }

            cellPoints->SetId(4, this->NumPoints + nAdditionalPoints);
            // decompose a face into quads in order (flipping the
            // decomposed face if owner)
            for (; vertI < nFaceJPoints; vertI += 2)
              {
              cellPoints->SetId(1, faceJPoints[(vertI + 1) % nFaceJPoints]);
              cellPoints->SetId(2, faceJPoints[vertI]);
              cellPoints->SetId(3, faceJPoints[vertI - 1]);

              // if the decomposed cell is the first one insert it to
              // the original position; or append to the decomposed cell
              // list otherwise
              if (insertDecomposedCell)
                {
                internalMesh->InsertNextCell(VTK_PYRAMID, 5,
                    cellPoints->GetPointer(0));
                insertDecomposedCell = false;
                }
              else
                {
                nAdditionalCells++;
                additionalCells->InsertNextTupleValue(
                    cellPoints->GetPointer(0));
                }
              }
            }
          else
            {
            const int flipNeighbor = (this->FaceOwner->GetValue(cellFacesJ)
                == cellId ? -1 : 1);

            cellPoints->SetId(0,
                faceJPoints[vertI == 2 ? 0 : nFaceJPoints - 1]);
            cellPoints->SetId(4, this->NumPoints + nAdditionalPoints);

            // decompose a face into quads in order (flipping the
            // decomposed face if owner)
            const int nQuadVerts = nFaceJPoints - 1 - nTris;
            for (; vertI < nQuadVerts; vertI += 2)
              {
              cellPoints->SetId(1, faceJPoints[vertI - flipNeighbor]);
              cellPoints->SetId(2, faceJPoints[vertI]);
              cellPoints->SetId(3, faceJPoints[vertI + flipNeighbor]);

              // if the decomposed cell is the first one insert it to
              // the original position; or append to the decomposed cell
              // list otherwise
              if (insertDecomposedCell)
                {
                internalMesh->InsertNextCell(VTK_PYRAMID, 5,
                    cellPoints->GetPointer(0));
                insertDecomposedCell = false;
                }
              else
                {
                nAdditionalCells++;
                additionalCells->InsertNextTupleValue(
                    cellPoints->GetPointer(0));
                }
              }

            // if the number of vertices is odd there's a triangle
            if (nTris)
              {
              if (flipNeighbor == -1)
                {
                cellPoints->SetId(1, faceJPoints[vertI]);
                cellPoints->SetId(2, faceJPoints[vertI - 1]);
                }
              else
                {
                cellPoints->SetId(1, faceJPoints[vertI - 1]);
                cellPoints->SetId(2, faceJPoints[vertI]);
                }
              cellPoints->SetId(3, this->NumPoints + nAdditionalPoints);

              if (insertDecomposedCell)
                {
                internalMesh->InsertNextCell(VTK_TETRA, 4,
                    cellPoints->GetPointer(0));
                insertDecomposedCell = false;
                }
              else
                {
                // set the 5th vertex number to -1 to distinguish a tetra cell
                cellPoints->SetId(4, -1);
                nAdditionalCells++;
                additionalCells->InsertNextTupleValue(
                    cellPoints->GetPointer(0));
                }
              }
            }
          }
        nAdditionalPoints++;
        this->AdditionalCellIds->InsertNextValue(cellId);
        this->NumAdditionalCells->InsertNextValue(nAdditionalCells);
        this->NumTotalAdditionalCells += nAdditionalCells;
        }
      else // don't decompose; use either VTK_POLYHEDRON or VTK_CONVEX_PONIT_SET
        {
        // get first face
        const int cellFaces0 = cellFaces[0];
        const int *baseFacePoints = facePoints[cellFaces0];
        const int nBaseFacePoints = facePoints.GetSize(cellFaces0);
#if VTK_FOAMFILE_USE_VTK_POLYHEDRON
        int nPoints = nBaseFacePoints, nPolyPoints = nBaseFacePoints + 1;
        if (nPoints > maxNPoints || nPolyPoints > maxNPolyPoints)
          {
          vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          cellPoints->Delete();
          polyPoints->Delete();
          return;
          }
        polyPoints->SetId(0, nBaseFacePoints);
        if (this->FaceOwner->GetValue(cellFaces0) == cellId)
          {
          // add first face to cell points
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            const int pointJ = baseFacePoints[j];
            cellPoints->SetId(j, pointJ);
            polyPoints->SetId(j + 1, pointJ);
            }
          }
        else
          {
          // if it is a _neighbor_ face flip the points
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            const int pointJ = baseFacePoints[nBaseFacePoints - 1 - j];
            cellPoints->SetId(j, pointJ);
            polyPoints->SetId(j + 1, pointJ);
            }
          }

        // loop through faces and create a list of all points
        // j = 1 skip baseFace
        for (int j = 1; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.GetSize(cellFacesJ);
          if (nPolyPoints >= maxNPolyPoints)
            {
            vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
            this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
            cellPoints->Delete();
            polyPoints->Delete();
            return;
            }
          polyPoints->SetId(nPolyPoints++, nFaceJPoints);
          int pointI, delta; // must be signed
          if (this->FaceOwner->GetValue(cellFacesJ) == cellId)
            {
            pointI = 0;
            delta = 1;
            }
          else
            {
            // if it is a _neighbor_ face flip the points
            pointI = static_cast<int>(nFaceJPoints) - 1;
            delta = -1;
            }
          for (size_t k = 0; k < nFaceJPoints; k++, pointI += delta)
            {
            const int faceJPointK = faceJPoints[pointI];
            bool foundDup = false;
            for (int l = 0; l < nPoints; l++)
              {
              if (cellPoints->GetId(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if (!foundDup)
              {
              if (nPoints >= maxNPoints)
                {
                vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
                this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
                cellPoints->Delete();
                polyPoints->Delete();
                return;
                }
              cellPoints->SetId(nPoints++, faceJPointK);
              }
            if (nPolyPoints >= maxNPolyPoints)
                {
                vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
                this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
                cellPoints->Delete();
                polyPoints->Delete();
                return;
                }
            polyPoints->SetId(nPolyPoints++, faceJPointK);
            }
          }

        // create the poly cell and insert it into the mesh
        internalMesh->InsertNextCell(VTK_POLYHEDRON, nPoints,
            cellPoints->GetPointer(0), nCellFaces, polyPoints->GetPointer(0));
#else
        int nPoints = nBaseFacePoints;
        if (nPoints > maxNPoints)
          {
          vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
          this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
          cellPoints->Delete();
          return;
          }
        if (this->FaceOwner->GetValue(cellFaces0) == cellId)
          {
          // if it is an owner face flip the points
          // not sure if flipping is necessary but do it anyway
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
            }
          }
        else
          {
          // add first face to cell points
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[j]);
            }
          }

        // loop through faces and create a list of all points
        // j = 1 skip baseFace
        for (int j = 1; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.GetSize(cellFacesJ);
          for (size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for (int l = 0; l < nPoints; l++)
              {
              if (cellPoints->GetId(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if (!foundDup)
              {
              if (nPoints >= maxNPoints)
                {
                vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
                this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
                cellPoints->Delete();
                return;
                }
              cellPoints->SetId(nPoints++, faceJPointK);
              }
            }
          }

        // create the poly cell and insert it into the mesh
        internalMesh->InsertNextCell(VTK_CONVEX_POINT_SET, nPoints,
            cellPoints->GetPointer(0));
#endif
        }
      }
    }
  cellPoints->Delete();
#if VTK_FOAMFILE_USE_VTK_POLYHEDRON
  polyPoints->Delete();
#endif
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::SetBlockName(vtkMultiBlockDataSet *blocks,
    unsigned int blockI, const char *name)
{
  blocks->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), name);
}

//-----------------------------------------------------------------------------
// derive cell types and create the internal mesh
vtkUnstructuredGrid *vtkOFFReaderPrivate::MakeInternalMesh(
    const vtkFoamIntVectorVector *cellsFaces,
    const vtkFoamIntVectorVector *facesPoints, vtkFloatArray *pointArray)
{
  // Create Mesh
  vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::New();
  internalMesh->Allocate(this->NumCells);

  if (this->Parent->GetDecomposePolyhedra())
    {
    // for polyhedral decomposition
    this->AdditionalCellIds = vtkIntArray::New();
    this->NumAdditionalCells = vtkIntArray::New();
    this->AdditionalCellPoints = new vtkFoamIdListVector;

    vtkIdTypeArray *additionalCells = vtkIdTypeArray::New();
    additionalCells->SetNumberOfComponents(5); // accommodates tetra or pyramid

    this->InsertCellsToGrid(internalMesh, cellsFaces, facesPoints, pointArray,
        additionalCells, NULL);

    // for polyhedral decomposition
    pointArray->Squeeze();
    this->AdditionalCellIds->Squeeze();
    this->NumAdditionalCells->Squeeze();
    additionalCells->Squeeze();

    // insert decomposed cells into mesh
    const int nComponents = additionalCells->GetNumberOfComponents();
    const int nAdditionalCells = additionalCells->GetNumberOfTuples();
    for (int i = 0; i < nAdditionalCells; i++)
      {
      if (additionalCells->GetComponent(i, 4) == -1)
        {
        internalMesh->InsertNextCell(VTK_TETRA, 4,
            additionalCells->GetPointer(i * nComponents));
        }
      else
        {
        internalMesh->InsertNextCell(VTK_PYRAMID, 5,
            additionalCells->GetPointer(i * nComponents));
        }
      }
    internalMesh->Squeeze();
    additionalCells->Delete();
    }
  else
    {
    this->InsertCellsToGrid(internalMesh, cellsFaces, facesPoints, pointArray,
        NULL, NULL);
    }

  // set the internal mesh points
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  internalMesh->SetPoints(points);
  points->Delete();

  return internalMesh;
}

//-----------------------------------------------------------------------------
// create surface mesh
vtkPolyData* vtkOFFReaderPrivate::MakeSurfaceMesh(
    const vtkStdString &meshDir, const vtkFoamIntVectorVector *facesPoints,
    vtkFloatArray *pointArray)
{
  if (this->BoundaryDict.size())
    {
    this->NumInternalFaces = this->BoundaryDict[0].StartFace;
    }
  else
    {
    this->NumInternalFaces = facesPoints->GetNumberOfElements();
    }

  // Create Mesh
  vtkPolyData* surfaceMesh = vtkPolyData::New();
  surfaceMesh->Allocate(this->NumInternalFaces);

  // Count the max number of points per face in mesh
  int maxNFacePoints = 0;
  for (int j = 0; j < this->NumInternalFaces; j++)
    {
    const int nFacePoints = facesPoints->GetSize(j);
    if (nFacePoints > maxNFacePoints)
      {
      maxNFacePoints = nFacePoints;
      }
    }

  vtkIdList *facePointsVtkId = vtkIdList::New();
  facePointsVtkId->SetNumberOfIds(maxNFacePoints);

  // Add faces to surfaceMesh
  this->InsertFacesToGrid(surfaceMesh, facesPoints, 0, this->NumInternalFaces,
      NULL, facePointsVtkId, NULL, false);

  const int nProcBoundaries = this->BoundaryDict.ProcBoundaries.size();
  if (nProcBoundaries > 0)
    {
    vtkFoamIOobject io(this->CasePath,
        this->Parent->GetIsSinglePrecisionBinary() != 0);
    vtkStdString faceProcPath(meshDir + "faceProcAddressing");
    if (!(io.Open(faceProcPath) || io.Open(faceProcPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      this->Parent->SetErrorCode(vtkErrorCode::CannotOpenFileError);
      return NULL;
      }

    vtkFoamEntryValue faceProcDict(NULL);
    try
      {
      faceProcDict.ReadNonuniformList<vtkFoamToken::LABELLIST,
        vtkFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }
    io.Close();

    // first pass: count the number of processor faces to be added
    vtkIntArray *faceProc = static_cast<vtkIntArray *>(faceProcDict.Ptr());
    int nNonNegFaces = 0;
    for (int procBoundI = 0; procBoundI < nProcBoundaries; procBoundI++)
      {
      vtkFoamBoundaryEntry &beI
          = this->BoundaryDict[this->BoundaryDict.ProcBoundaries[procBoundI]];
      const int startFace = beI.StartFace, endFace = startFace + beI.NFaces;
      for (int faceI = startFace; faceI < endFace; faceI++)
        {
        if (faceProc->GetValue(faceI) >= 0)
          {
          nNonNegFaces++;
          const int nFacePoints = facesPoints->GetSize(faceI);
          if (nFacePoints > maxNFacePoints)
            {
            maxNFacePoints = nFacePoints;
            }
          }
        }
      }

    if (nNonNegFaces > 0)
      {
      // second pass: fill the processor face id vector and insert faces
      facePointsVtkId->SetNumberOfIds(maxNFacePoints);
      this->ProcessorFaces
          = new vtkFoamIntVectorVector(nProcBoundaries, nNonNegFaces);
      nNonNegFaces = 0;
      for (int procBoundI = 0; procBoundI < nProcBoundaries; procBoundI++)
        {
        int *procFaces = this->ProcessorFaces->SetIndex(procBoundI,
            nNonNegFaces);
        vtkFoamBoundaryEntry &beI
            = this->BoundaryDict[this->BoundaryDict.ProcBoundaries[procBoundI]];
        const int startFace = beI.StartFace, endFace = startFace + beI.NFaces;
        int nProcFaces = 0;
        for (int faceI = startFace; faceI < endFace; faceI++)
          {
          if (faceProc->GetValue(faceI) >= 0)
            {
            procFaces[nProcFaces++] = faceI;
            }
          }
        this->InsertFacesToGrid(surfaceMesh, facesPoints, 0, nProcFaces, NULL,
            facePointsVtkId, procFaces, false);
        for (int faceI = 0; faceI < nProcFaces; faceI++)
          {
          procFaces[faceI] -= startFace;
          }
        nNonNegFaces += nProcFaces;
        }
      this->ProcessorFaces->GetIndices()->SetValue(nProcBoundaries,
          nNonNegFaces);
      }
    faceProc->Delete();
    }

  facePointsVtkId->Delete();

  // Add marker for "Show Region Names" operation
  vtkIntArray *bt = vtkIntArray::New();
  bt->SetNumberOfTuples(1);
  bt->SetValue(0, vtkFoamBoundaryEntry::INTERNAL);
  bt->SetName("BoundaryType");
  surfaceMesh->GetFieldData()->AddArray(bt);
  bt->FastDelete();

  // Set the internal mesh points
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  surfaceMesh->SetPoints(points);
  points->Delete();

  return surfaceMesh;
}

//-----------------------------------------------------------------------------
// insert faces to grid
bool vtkOFFReaderPrivate::InsertFacesToGrid(vtkPolyData *boundaryMesh,
    const vtkFoamIntVectorVector *facesPoints, int startFace, int endFace,
    vtkIntArray *boundaryPointMap, vtkIdList *facePointsVtkId,
    const int *labels, bool isLookupValue)
{
  vtkPolyData &bm = *boundaryMesh;

  for (int j = startFace; j < endFace; j++)
    {
    int faceId;
    if (labels == NULL)
      {
      faceId = j;
      }
    else
      {
      faceId = labels[j];
      if (faceId >= facesPoints->GetNumberOfElements())
        {
        vtkErrorMacro(<<"faceLabels id " << faceId
            << " exceeds the number of faces "
            << facesPoints->GetNumberOfElements());
        this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
        return false;
        }
      }
    const int *facePoints = facesPoints->operator[](faceId);
    vtkIdType nFacePoints = facesPoints->GetSize(faceId);

    if (isLookupValue)
      {
      for (vtkIdType k = 0; k < nFacePoints; k++)
        {
        facePointsVtkId->SetId(k, boundaryPointMap->LookupValue(facePoints[k]));
        }
      }
    else
      {
      if (boundaryPointMap)
        {
        for (vtkIdType k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, boundaryPointMap->GetValue(facePoints[k]));
          }
        }
      else
        {
        for (vtkIdType k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, facePoints[k]);
          }
        }
      }

    // triangle
    if (nFacePoints == 3)
      {
      bm.InsertNextCell(VTK_TRIANGLE, 3, facePointsVtkId->GetPointer(0));
      }
    // quad
    else if (nFacePoints == 4)
      {
      bm.InsertNextCell(VTK_QUAD, 4, facePointsVtkId->GetPointer(0));
      }
    // polygon
    else
      {
      bm.InsertNextCell(VTK_POLYGON, nFacePoints,
          facePointsVtkId->GetPointer(0));
      }
    }
  return true;
}

//-----------------------------------------------------------------------------
// returns requested boundary meshes
vtkMultiBlockDataSet *vtkOFFReaderPrivate::MakeBoundaryMesh(
    const vtkFoamIntVectorVector *facesPoints, vtkFloatArray* pointArray)
{
  const vtkIdType nBoundaries = static_cast<vtkIdType>(this->BoundaryDict.size());

  // do a consistency check of BoundaryDict
  int previousEndFace = -1;
  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int startFace = beI.StartFace;
    const int nFaces = beI.NFaces;
    if (startFace < 0 || nFaces < 0)
      {
      vtkErrorMacro(<<"Neither of startFace " << startFace << " nor nFaces "
          << nFaces << " can be nagative for patch " << beI.BoundaryName.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }
    if (previousEndFace >= 0 && previousEndFace != startFace)
      {
      vtkErrorMacro(<<"The end face number " << previousEndFace - 1
          << " of patch "
          << this->BoundaryDict[boundaryI - 1].BoundaryName.c_str()
          << " is not consistent with the start face number " << startFace
          << " of patch " << beI.BoundaryName.c_str());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return NULL;
      }
    previousEndFace = startFace + nFaces;
    }
  if (previousEndFace > facesPoints->GetNumberOfElements())
    {
    vtkErrorMacro(<<"The end face number " << previousEndFace - 1
        << " of the last patch "
        << this->BoundaryDict[nBoundaries - 1].BoundaryName.c_str()
        << " exceeds the number of faces " << facesPoints->GetNumberOfElements());
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return NULL;
    }

  vtkMultiBlockDataSet *boundaryMesh = vtkMultiBlockDataSet::New();

  if (this->Parent->GetCreateCellToPoint())
    {
    const int boundaryStartFace =
        (this->BoundaryDict.size() > 0 ? this->BoundaryDict[0].StartFace : 0);
    this->AllBoundaries = vtkPolyData::New();
    this->AllBoundaries->Allocate(facesPoints->GetNumberOfElements()
        - boundaryStartFace);
    }
  this->BoundaryPointMap = new vtkFoamIntArrayVector;

  vtkIntArray *nBoundaryPointsList = vtkIntArray::New();
  nBoundaryPointsList->SetNumberOfValues(nBoundaries);

  // count the max number of points per face and the number of points
  // (with duplicates) in mesh
  int maxNFacePoints = 0;
  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const int startFace = this->BoundaryDict[boundaryI].StartFace;
    const int endFace = startFace + this->BoundaryDict[boundaryI].NFaces;
    int nPoints = 0;
    for (int j = startFace; j < endFace; j++)
      {
      const int nFacePoints = facesPoints->GetSize(j);
      nPoints += nFacePoints;
      if (nFacePoints > maxNFacePoints)
        {
        maxNFacePoints = nFacePoints;
        }
      }
    nBoundaryPointsList->SetValue(boundaryI, nPoints);
    }

  // aloocate array for converting int vector to vtkIdType List:
  // workaround for 64bit machines
  vtkIdList *facePointsVtkId = vtkIdList::New();
  facePointsVtkId->SetNumberOfIds(maxNFacePoints);

  // create initial internal point list: set all points to -1
  if (this->Parent->GetCreateCellToPoint())
    {
    this->InternalPoints = vtkIntArray::New();
    this->InternalPoints->SetNumberOfValues(this->NumPoints);
    for (int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      this->InternalPoints->SetValue(pointI, -1);
      }

    // mark boundary points as 0
    for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
      {
      const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
      if (beI.BoundaryType == vtkFoamBoundaryEntry::PHYSICAL
          || beI.BoundaryType == vtkFoamBoundaryEntry::PROCESSOR)
        {
        const int startFace = beI.StartFace;
        const int endFace = startFace + beI.NFaces;

        for (int j = startFace; j < endFace; j++)
          {
          const int *facePoints = facesPoints->operator[](j);
          const int nFacePoints = facesPoints->GetSize(j);
          for (int k = 0; k < nFacePoints; k++)
            {
            this->InternalPoints->SetValue(facePoints[k], 0);
            }
          }
        }
      }
    }

  int nAllBoundaryPoints = 0;
  vtkstd::vector<vtkstd::vector<int> > procCellList;
  vtkIntArray *pointTypes = NULL;

  if (this->Parent->GetCreateCellToPoint())
    {
    // create global to AllBoundaries point map
    for (int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      if (this->InternalPoints->GetValue(pointI) == 0)
        {
        this->InternalPoints->SetValue(pointI, nAllBoundaryPoints);
        nAllBoundaryPoints++;
        }
      }

    if (this->ProcessorName != "")
      {
      // initialize physical-processor boundary shared point list
      procCellList.resize(nAllBoundaryPoints);
      pointTypes = vtkIntArray::New();
      pointTypes->SetNumberOfTuples(nAllBoundaryPoints);
      for (int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        pointTypes->SetValue(pointI, 0);
        }
      }
    }

  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int nFaces = beI.NFaces;
    const int startFace = beI.StartFace;
    const int endFace = startFace + nFaces;

    if (this->Parent->GetCreateCellToPoint() && (beI.BoundaryType
        == vtkFoamBoundaryEntry::PHYSICAL || beI.BoundaryType
        == vtkFoamBoundaryEntry::PROCESSOR))
      {
      // add faces to AllBoundaries
      this->InsertFacesToGrid(this->AllBoundaries, facesPoints, startFace,
          endFace, this->InternalPoints, facePointsVtkId, NULL, false);

      if (this->ProcessorName != "")
        {
        // mark belonging boundary types and, if PROCESSOR, cell numbers
        const int abStartFace = beI.AllBoundariesStartFace;
        const int abEndFace = abStartFace + beI.NFaces;
        for (int faceI = abStartFace; faceI < abEndFace; faceI++)
          {
          vtkIdType nPoints;
          vtkIdType *points;
          this->AllBoundaries->GetCellPoints(faceI, nPoints, points);
          if (beI.BoundaryType == vtkFoamBoundaryEntry::PHYSICAL)
            {
            for (int pointI = 0; pointI < nPoints; pointI++)
              {
              *pointTypes->GetPointer(points[pointI])
                  |= vtkFoamBoundaryEntry::PHYSICAL;
              }
            }
          else // PROCESSOR
            {
            for (int pointI = 0; pointI < nPoints; pointI++)
              {
              const int pointJ = points[pointI];
              *pointTypes->GetPointer(pointJ)
                  |= vtkFoamBoundaryEntry::PROCESSOR;
              procCellList[pointJ].push_back(faceI);
              }
            }
          }
        }
      }

    // skip below if inactive
    if (!beI.IsActive)
      {
      continue;
      }

    // create the mesh
    const unsigned int activeBoundaryI = boundaryMesh->GetNumberOfBlocks();
    vtkPolyData *bm = vtkPolyData::New();
    boundaryMesh->SetBlock(activeBoundaryI, bm);

    // set the name of boundary
    this->SetBlockName(boundaryMesh, activeBoundaryI, beI.BoundaryName.c_str());

    bm->Allocate(nFaces);
    const int nBoundaryPoints = nBoundaryPointsList->GetValue(boundaryI);

    // create global to boundary-local point map and boundary points
    vtkIntArray *boundaryPointList = vtkIntArray::New();
    boundaryPointList->SetNumberOfValues(nBoundaryPoints);
    int pointI = 0;
    for (int j = startFace; j < endFace; j++)
      {
      const int *facePoints = facesPoints->operator[](j);
      int nFacePoints = facesPoints->GetSize(j);
      for (int k = 0; k < nFacePoints; k++)
        {
        boundaryPointList->SetValue(pointI, facePoints[k]);
        pointI++;
        }
      }
    vtkSortDataArray::Sort(boundaryPointList);
    this->BoundaryPointMap->push_back(vtkIntArray::New());
    vtkIntArray& bpMap = *this->BoundaryPointMap->back();
    vtkFloatArray *boundaryPointArray = vtkFloatArray::New();
    boundaryPointArray->SetNumberOfComponents(3);
    int oldPointJ = -1;
    for (int j = 0; j < nBoundaryPoints; j++)
      {
      const int pointJ = boundaryPointList->GetValue(j);
      if (pointJ != oldPointJ)
        {
        oldPointJ = pointJ;
        boundaryPointArray->InsertNextTuple(pointArray->GetPointer(3 * pointJ));
        bpMap.InsertNextValue(pointJ);
        }
      }
    boundaryPointArray->Squeeze();
    bpMap.Squeeze();
    boundaryPointList->Delete();
    vtkPoints *boundaryPoints = vtkPoints::New();
    boundaryPoints->SetData(boundaryPointArray);
    boundaryPointArray->Delete();

    // set points for boundary
    bm->SetPoints(boundaryPoints);
    boundaryPoints->Delete();

    // insert faces to boundary mesh
    this->InsertFacesToGrid(bm, facesPoints, startFace, endFace, &bpMap,
        facePointsVtkId, NULL, true);

    vtkIntArray *bt = vtkIntArray::New();
    bt->SetNumberOfTuples(1);
    bt->SetValue(0, beI.BoundaryType);
    bt->SetName("BoundaryType");
    bm->GetFieldData()->AddArray(bt);
    bt->FastDelete();
    if (beI.BoundaryType == vtkFoamBoundaryEntry::PROCESSOR)
      {
      vtkIntArray *procNo = vtkIntArray::New();
      procNo->SetNumberOfTuples(1);
      procNo->SetValue(0, beI.MyProcNo);
      procNo->SetName("myProcNo");
      bm->GetFieldData()->AddArray(procNo);
      procNo->FastDelete();

      procNo = vtkIntArray::New();
      procNo->SetNumberOfTuples(1);
      procNo->SetValue(0, beI.NeighbProcNo);
      procNo->SetName("neighbProcNo");
      bm->GetFieldData()->AddArray(procNo);
      procNo->FastDelete();
      }

    bm->Delete();
    bpMap.ClearLookup();
    }

  nBoundaryPointsList->Delete();
  facePointsVtkId->Delete();

  if (this->Parent->GetCreateCellToPoint())
    {
    this->AllBoundaries->Squeeze();
    this->AllBoundariesPointMap = vtkIntArray::New();
    vtkIntArray &abpMap = *this->AllBoundariesPointMap;
    abpMap.SetNumberOfValues(nAllBoundaryPoints);

    // create lists of internal points and AllBoundaries points
    int nInternalPoints = 0;
    for (int pointI = 0, allBoundaryPointI = 0; pointI < this->NumPoints; pointI++)
      {
      const int globalPointId = this->InternalPoints->GetValue(pointI);
      if (globalPointId == -1)
        {
        this->InternalPoints->SetValue(nInternalPoints, pointI);
        nInternalPoints++;
        }
      else
        {
        abpMap.SetValue(allBoundaryPointI, pointI);
        allBoundaryPointI++;
        }
      }
    // shrink to the number of internal points
    if (nInternalPoints > 0)
      {
      this->InternalPoints->Resize(nInternalPoints);
      }
    else
      {
      this->InternalPoints->Delete();
      this->InternalPoints = NULL;
      }

    // set dummy vtkPoints to tell the grid the number of points
    // (otherwise GetPointCells will crash)
    vtkPoints *allBoundaryPoints = vtkPoints::New();
    allBoundaryPoints->SetNumberOfPoints(abpMap.GetNumberOfTuples());
    this->AllBoundaries->SetPoints(allBoundaryPoints);
    allBoundaryPoints->Delete();

    if (this->ProcessorName != "")
      {
      // remove links to processor boundary faces from point-to-cell
      // links of physical-processor shared points to avoid cracky seams
      // on fixedValue-type boundaries which are noticeable when all the
      // decomposed meshes are appended
      this->AllBoundaries->BuildLinks();
      for (int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        if (pointTypes->GetValue(pointI) == (vtkFoamBoundaryEntry::PHYSICAL
            | vtkFoamBoundaryEntry::PROCESSOR))
          {
          const vtkstd::vector<int> &procCells = procCellList[pointI];
          for (size_t cellI = 0; cellI < procCellList[pointI].size(); cellI++)
            {
            this->AllBoundaries->RemoveReferenceToCell(pointI, procCells[cellI]);
            }
          // omit reclaiming memory as the possibly recovered size should
          // not typically be so large
          }
        }
      pointTypes->Delete();
      }
    }

  return boundaryMesh;
}

//-----------------------------------------------------------------------------
// truncate face owner to have only boundary face info
void vtkOFFReaderPrivate::TruncateFaceOwner()
{
  const int boundaryStartFace =
      this->BoundaryDict.size() > 0 ? this->BoundaryDict[0].StartFace
          : this->FaceOwner->GetNumberOfTuples();
  // all the boundary faces
  const int nBoundaryFaces = this->FaceOwner->GetNumberOfTuples()
      - boundaryStartFace;
  memmove(this->FaceOwner->GetPointer(0),
      this->FaceOwner->GetPointer(boundaryStartFace), sizeof(int)
          * nBoundaryFaces);
  this->FaceOwner->Resize(nBoundaryFaces);
}

//-----------------------------------------------------------------------------
// this is necessary due to the strange vtkDataArrayTemplate::Resize()
// implementation when the array size is to be extended
template <typename T1, typename T2> bool vtkOFFReaderPrivate::ExtendArray(
    T1 *array, const int nTuples)
{
  const int newSize = nTuples * array->GetNumberOfComponents();
  void *ptr = malloc(newSize * array->GetDataTypeSize());
  if (ptr == NULL)
    {
    return false;
    }
  memmove(ptr, array->GetVoidPointer(0), array->GetDataSize()
      * array->GetDataTypeSize());
  array->SetArray(static_cast<T2 *>(ptr), newSize, 0);
  return true;
}

#if 0
//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::CalculateReciprocalDelta(
  const vtkFoamIntVectorVector *facesPoints)
{
  const int nBoundaries = this->BoundaryDict.size();
  vtkFloatArray *allPoints
    = vtkFloatArray::SafeDownCast(this->InternalMesh->GetPoints()->GetData());

  this->ReciprocalDelta = new vtkFoamFloatArrayVector;
  for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    if(beI.BoundaryType != vtkFoamBoundaryEntry::PHYSICAL)
      {
      this->ReciprocalDelta->push_back(NULL);
      continue;
      }

    this->ReciprocalDelta->push_back(vtkFloatArray::New());
    vtkFloatArray *rd = this->ReciprocalDelta->back();
    rd->SetNumberOfTuples(beI.NFaces);
    const int boundaryStartFace = this->BoundaryDict[0].StartFace;
    const int startFace = beI.StartFace;
    const int endFace = startFace + beI.NFaces;
    for(int faceI = startFace; faceI < endFace; faceI++)
      {
      // calculate patchInternal cell centroid
      const int cellId = this->FaceOwner->GetValue(faceI - boundaryStartFace);
      const int polyCellId = this->AdditionalCellIds->LookupValue(cellId);
      float cellCentroid0, cellCentroid1, cellCentroid2;
      if(this->Parent->GetDecomposePolyhedra()
        && this->AdditionalCellPoints != NULL && polyCellId >= 0)
        {
        const float *cellCentroidPtr
          = allPoints->GetPointer(3 * (this->NumPoints + polyCellId));
        cellCentroid0 = cellCentroidPtr[0];
        cellCentroid1 = cellCentroidPtr[1];
        cellCentroid2 = cellCentroidPtr[2];
        }
      else
        {
        cellCentroid0 = cellCentroid1 = cellCentroid2 = 0.0F;
        vtkIdType nPoints, *pointIds;
        this->InternalMesh->GetCellPoints(cellId, nPoints, pointIds);
        for(int pointI = 0; pointI < nPoints; pointI++)
          {
          const float *point = allPoints->GetPointer(3 * pointIds[pointI]);
          cellCentroid0 += point[0];
          cellCentroid1 += point[1];
          cellCentroid2 += point[2];
          }
        const float weight = 1.0F / static_cast<float>(nPoints);
        cellCentroid0 *= weight;
        cellCentroid1 *= weight;
        cellCentroid2 *= weight;
        }

      // calculate face centroid
      float faceCentroid0 = 0.0F, faceCentroid1 = 0.0F, faceCentroid2 = 0.0F;
      const int *facePoints = facesPoints->operator[](faceI);
      const int nFacePoints = facesPoints->GetSize(faceI);
      for(int pointI = 0; pointI < nFacePoints; pointI++)
        {
        const float *point = allPoints->GetPointer(3 * facePoints[pointI]);
        faceCentroid0 += point[0];
        faceCentroid1 += point[1];
        faceCentroid2 += point[2];
        }
      const float weight = 1.0F / static_cast<float>(nFacePoints);
      faceCentroid0 *= weight;
      faceCentroid1 *= weight;
      faceCentroid2 *= weight;

      // calculate face normal
      const float *point0 = allPoints->GetPointer(3 * facePoints[0]);
      const float *point1 = allPoints->GetPointer(3 * facePoints[1]);
      const float v00 = point0[0] - point1[0];
      const float v01 = point0[1] - point1[1];
      const float v02 = point0[2] - point1[2];
      const float *point2 = allPoints->GetPointer(3 * facePoints[2]);
      const float v10 = point2[0] - point1[0];
      const float v11 = point2[1] - point1[1];
      const float v12 = point2[2] - point1[2];
      float n0 = v11 * v02 - v12 * v01;
      float n1 = v12 * v00 - v10 * v02;
      float n2 = v10 * v01 - v11 * v00;
      float nWeight = sqrt(n0 * n0 + n1 * n1 + n2 * n2);

      const float delta0 = faceCentroid0 - cellCentroid0;
      const float delta1 = faceCentroid1 - cellCentroid1;
      const float delta2 = faceCentroid2 - cellCentroid2;
      float delta;
      if(nWeight == 0.0F)
        {
        // use raw distance
        delta = sqrt(delta0 * delta0 + delta1 * delta1 + delta2 * delta2);
        }
      else
        {
        // use inner product
        delta = (delta0 * n0 + delta1 * n1 + delta2 * n2) / nWeight;
        }
      rd->SetValue(faceI - beI.StartFace, delta == 0.0F ? delta : 1.0F / delta);
      }
    }
  this->AdditionalCellIds->ClearLookup();
}
#endif

//-----------------------------------------------------------------------------
// move polyhedral cell centroids
vtkPoints *vtkOFFReaderPrivate::MoveInternalMesh(
    vtkUnstructuredGrid *internalMesh, vtkFloatArray *pointArray)
{
  if (this->Parent->GetDecomposePolyhedra())
    {
    const vtkIdType nAdditionalCells = static_cast<vtkIdType>(this->AdditionalCellPoints->size());
    this->ExtendArray<vtkFloatArray, float>(pointArray, this->NumPoints
        + nAdditionalCells);
    for (int i = 0; i < nAdditionalCells; i++)
      {
      vtkIdList *polyCellPoints = this->AdditionalCellPoints->operator[](i);
      float centroid[3];
      centroid[0] = centroid[1] = centroid[2] = 0.0F;
      const vtkIdType nCellPoints = polyCellPoints->GetNumberOfIds();
      for (vtkIdType j = 0; j < nCellPoints; j++)
        {
        float *pointK = pointArray->GetPointer(3 * polyCellPoints->GetId(j));
        centroid[0] += pointK[0];
        centroid[1] += pointK[1];
        centroid[2] += pointK[2];
        }
      const float weight = (nCellPoints ? 1.0F
          / static_cast<float>(nCellPoints) : 0.0F);
      centroid[0] *= weight;
      centroid[1] *= weight;
      centroid[2] *= weight;
      pointArray->InsertTuple(this->NumPoints + i, centroid);
      }
    }
  if (internalMesh->GetPoints()->GetNumberOfPoints() != pointArray->GetNumberOfTuples())
    {
    vtkErrorMacro(<< "The numbers of points for old points "
        << internalMesh->GetPoints()->GetNumberOfPoints() << " and new points"
        << pointArray->GetNumberOfTuples() << " don't match");
    this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
    return NULL;
    }

  // instantiate the points class
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  internalMesh->SetPoints(points);
  return points;
}

//-----------------------------------------------------------------------------
// move internal surface mesh
int vtkOFFReaderPrivate::MoveSurfaceMesh(vtkPolyData *surfaceMesh,
    vtkFloatArray *pointArray)
{
  if (surfaceMesh->GetPoints()->GetNumberOfPoints()
      != pointArray->GetNumberOfTuples())
    {
    vtkErrorMacro(<< "The numbers of points for old points "
        << surfaceMesh->GetPoints()->GetNumberOfPoints() << " and new points"
        << pointArray->GetNumberOfTuples() << " don't match");

    return 0;
    }

  // instantiate the points class
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  surfaceMesh->SetPoints(points);
  points->Delete();
  return 1;
}

//-----------------------------------------------------------------------------
// move boundary points
void vtkOFFReaderPrivate::MoveBoundaryMesh(
    vtkMultiBlockDataSet *boundaryMesh, vtkFloatArray *pointArray)
{
  for (vtkIdType boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<vtkIdType>(this->BoundaryDict.size()); boundaryI++)
    {
    if (this->BoundaryDict[boundaryI].IsActive)
      {
      vtkIntArray *bpMap = this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nBoundaryPoints = bpMap->GetNumberOfTuples();
      vtkFloatArray *boundaryPointArray = vtkFloatArray::New();
      boundaryPointArray->SetNumberOfComponents(3);
      boundaryPointArray->SetNumberOfTuples(nBoundaryPoints);
      for (int pointI = 0; pointI < nBoundaryPoints; pointI++)
        {
        boundaryPointArray->SetTuple(pointI, bpMap->GetValue(pointI),
            pointArray);
        }
      vtkPoints *boundaryPoints = vtkPoints::New();
      boundaryPoints->SetData(boundaryPointArray);
      boundaryPointArray->Delete();
      vtkPolyData::SafeDownCast(boundaryMesh->GetBlock(activeBoundaryI))
      ->SetPoints(boundaryPoints);
      boundaryPoints->Delete();
      activeBoundaryI++;
      }
    }
}

//-----------------------------------------------------------------------------
// as of now the function does not do interpolation, but do just averaging.
void vtkOFFReaderPrivate::InterpolateCellToPoint(vtkFloatArray *pData,
    vtkFloatArray *iData, vtkPointSet *mesh, vtkIntArray *pointList,
    const int nPoints)
{
  if (nPoints == 0)
    {
    return;
    }

  // a dummy call to let GetPointCells() build the cell links if still not built
  // (not using BuildLinks() since it always rebuild links)
  vtkIdList *pointCells = vtkIdList::New();
  mesh->GetPointCells(0, pointCells);
  pointCells->Delete();

  // since vtkPolyData and vtkUnstructuredGrid do not share common
  // overloaded GetCellLink() or GetPointCells() functions we have to
  // do a tedious task
  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(mesh);
  vtkPolyData *pd = vtkPolyData::SafeDownCast(mesh);
  vtkCellLinks *cl = NULL;
  if (ug)
    {
    cl = ug->GetCellLinks();
    }

  const int nComponents = iData->GetNumberOfComponents();

  if (nComponents == 1)
    {
    // a special case with the innermost componentI loop unrolled
    float *tuples = iData->GetPointer(0);
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variable for precision
      double interpolatedValue = 0.0;
      for (int cellI = 0; cellI < nCells; cellI++)
        {
        interpolatedValue += tuples[cells[cellI]];
        }
      interpolatedValue = (nCells ? interpolatedValue
          / static_cast<double>(nCells) : 0.0);
      pData->SetValue(pI, interpolatedValue);
      }
    }
  else if (nComponents == 3)
    {
    // a special case with the innermost componentI loop unrolled
    float *pDataPtr = pData->GetPointer(0);
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      double summedValue0 = 0.0, summedValue1 = 0.0, summedValue2 = 0.0;

      // hand unrolling
      for (int cellI = 0; cellI < nCells; cellI++)
        {
        const float *tuple = iData->GetPointer(3 * cells[cellI]);
        summedValue0 += tuple[0];
        summedValue1 += tuple[1];
        summedValue2 += tuple[2];
        }

      float *interpolatedValue = &pDataPtr[3 * pI];
      interpolatedValue[0] = weight * summedValue0;
      interpolatedValue[1] = weight * summedValue1;
      interpolatedValue[2] = weight * summedValue2;
      }
    }
  else
    {
    float *pDataPtr = pData->GetPointer(0);
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      float *interpolatedValue = &pDataPtr[nComponents * pI];
      // a bit strange loop order but this works fastest
      for (int componentI = 0; componentI < nComponents; componentI++)
        {
        const float *tuple = iData->GetPointer(componentI);
        double summedValue = 0.0;
        for (int cellI = 0; cellI < nCells; cellI++)
          {
          summedValue += tuple[nComponents * cells[cellI]];
          }
        interpolatedValue[componentI] = weight * summedValue;
        }
      }
    }
}

//-----------------------------------------------------------------------------
bool vtkOFFReaderPrivate::ReadFieldFile(vtkFoamIOobject *ioPtr,
    vtkFoamDict *dictPtr, const vtkStdString &varName,
    vtkDataArraySelection *selection)
{
  const vtkStdString varPath(this->CurrentTimeRegionPath() + "/" + varName);

  // open the file
  vtkFoamIOobject &io = *ioPtr;
  if (!io.Open(varPath))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str());
    return false;
    }

  // if the variable is disabled on selection panel then skip it
  if (selection->ArrayExists(io.GetObjectName().c_str()) && !selection->ArrayIsEnabled(io.GetObjectName().c_str()))
    {
    return false;
    }

  // read the field file into dictionary
  vtkFoamDict &dict = *dictPtr;
  if (!dict.Read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
    return false;
    }

  if (dict.GetType() != vtkFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"File " << io.GetFileName().c_str()
        << "is not valid as a field file");
    return false;
    }
  return true;
}

//-----------------------------------------------------------------------------
vtkFloatArray *vtkOFFReaderPrivate::FillField(vtkFoamEntry *entryPtr,
    vtkFoamEntry *refEntryPtr, int nElements, vtkFoamIOobject *ioPtr,
    const vtkStdString &fieldType, const bool isUniformFixedValueBC)
    // the isUniformFixedValueBC argument can be determined in compile-time
    // thus better use a template, which unfortunately is not possible due to
    // broken gcc-3.3
{
  vtkFloatArray *data;
  vtkFoamEntry &entry = *entryPtr;
  const vtkStdString &className = ioPtr->GetClassName();

  // refEntryType and nRefComponents have to be initialized to prevent warning
  vtkFoamToken::tokenType refEntryType = vtkFoamToken::UNDEFINED;
  int nRefComponents = 0;
  float refTuple[9];
  if (refEntryPtr != NULL)
    {
    if (refEntryPtr->FirstValue().GetType() == vtkFoamToken::SCALAR
        || refEntryPtr->FirstValue().GetType() == vtkFoamToken::LABEL)
      {
      refEntryType = vtkFoamToken::SCALAR;
      nRefComponents = 1;
      refTuple[0] = refEntryPtr->ToFloat();
      }
    else if (refEntryPtr->FirstValue().GetType() == vtkFoamToken::LABELLIST)
      {
      refEntryType = vtkFoamToken::SCALARLIST;
      vtkIntArray &ll = refEntryPtr->LabelList();
      nRefComponents = ll.GetNumberOfTuples();
      for (int componentI = 0; componentI < nRefComponents; componentI++)
        {
        refTuple[componentI] = static_cast<float>(ll.GetValue(componentI));
        }
      }
    else if (refEntryPtr->FirstValue().GetType() == vtkFoamToken::SCALARLIST)
      {
      refEntryType = vtkFoamToken::SCALARLIST;
      vtkFloatArray &sl = refEntryPtr->ScalarList();
      nRefComponents = sl.GetNumberOfTuples();
      for (int componentI = 0; componentI < nRefComponents; componentI++)
        {
        refTuple[componentI] = sl.GetValue(componentI);
        }
      }
    else
      {
      vtkErrorMacro(<< "Wrong referenceLevel type");
      return NULL;
      }
    }

  // the isUniformFixedValueBC argument is for uniformFixedValue B.C.
  if (entry.FirstValue().GetIsUniform() == vtkFoamEntryValue::UNIFORM
      || isUniformFixedValueBC)
    {
    if (entry.FirstValue().GetType() == vtkFoamToken::SCALAR || entry.FirstValue().GetType() == vtkFoamToken::LABEL)
      {
      float num = entry.ToFloat();
      if (refEntryPtr != NULL)
        {
        if (refEntryType == vtkFoamToken::SCALAR)
          {
          num += refTuple[0];
          }
        else
          {
          vtkErrorMacro(<<"Wrong referenceLevel type for uniform scalar field");
          return NULL;
          }
        }
      data = vtkFloatArray::New();
      data->SetNumberOfValues(nElements);
      for (int i = 0; i < nElements; i++)
        {
        data->SetValue(i, num);
        }
      }
    else
      {
      float tuple[9];
      int nComponents;
      // have to determine the type of vector
      if (entry.FirstValue().GetType() == vtkFoamToken::LABELLIST)
        {
        vtkIntArray &ll = entry.LabelList();
        nComponents = ll.GetNumberOfTuples();
        for (int componentI = 0; componentI < nComponents; componentI++)
          {
          tuple[componentI] = static_cast<float>(ll.GetValue(componentI));
          }
        }
      else if (entry.FirstValue().GetType() == vtkFoamToken::SCALARLIST)
        {
        vtkFloatArray& sl = entry.ScalarList();
        nComponents = sl.GetSize();
        for (int componentI = 0; componentI < nComponents; componentI++)
          {
          tuple[componentI] = sl.GetValue(componentI);
          }
        }
      else
        {
        vtkErrorMacro(<< "Wrong list type for uniform field");
        return NULL;
        }

      if (refEntryPtr != NULL)
        {
        if (refEntryType == vtkFoamToken::SCALARLIST && nRefComponents == nComponents)
          {
          for (int componentI = 0; componentI < nComponents; componentI++)
            {
            tuple[componentI] += refTuple[componentI];
            }
          }
        else
          {
          vtkErrorMacro(<< "Wrong referenceLevel type");
          return NULL;
          }
        }

      if ((fieldType == "SphericalTensorField" && nComponents == 1)
          || (fieldType == "VectorField" && nComponents == 3) || (fieldType
          == "SymmTensorField" && nComponents == 6) || (fieldType
          == "TensorField" && nComponents == 9))
        {
        data = vtkFloatArray::New();
        data->SetNumberOfComponents(nComponents);
        data->SetNumberOfTuples(nElements);
#if VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS
        // swap the components of symmTensor to match the component
        // names in paraview
        // OpenFOAM: XX, XY, XZ, YY, YZ, ZZ
        // ParaView: XX, YY, ZZ, XY, YZ, XZ
        // if VTK is sufficiently new then use SetComponentName() instead
        if (nComponents == 6)
          {
          const float symxy = tuple[1], symxz = tuple[2], symyy = tuple[3];
          const float symzz = tuple[5];
          tuple[1] = symyy;
          tuple[2] = symzz;
          tuple[3] = symxy;
          tuple[5] = symxz;
          }
#endif
        for (int i = 0; i < nElements; i++)
          {
          data->SetTuple(i, tuple);
          }
        }
      else
        {
        vtkErrorMacro(<< "Number of components and field class doesn't match "
                      << "for " << ioPtr->GetFileName().c_str() << ". class = " << className.c_str()
                      << ", nComponents = " << nComponents);
        return NULL;
        }
      }
    }
  else if (entry.FirstValue().GetIsUniform() == vtkFoamEntryValue::NONUNIFORM)
    {
    if ((fieldType == "ScalarField" && entry.FirstValue().GetType() == vtkFoamToken::SCALARLIST) || ((fieldType
        == "VectorField" || fieldType == "SphericalTensorField" || fieldType
        == "SymmTensorField" || fieldType == "TensorField")
        && entry.FirstValue().GetType() == vtkFoamToken::VECTORLIST))
      {
      const int nTuples = entry.ScalarList().GetNumberOfTuples();
      if (nTuples != nElements)
        {
        vtkErrorMacro(<<"Number of cells/faces/points in mesh and field don't "
            << "match: mesh = " << nElements << ", field = " << nTuples);
        return NULL;
        }
      data = static_cast<vtkFloatArray *>(entry.Ptr());

      if (refEntryPtr != NULL)
        {
        const int nComponents = data->GetNumberOfComponents();
        if ((refEntryType == vtkFoamToken::SCALAR
            && entry.FirstValue().GetType() == vtkFoamToken::SCALARLIST)
            || (refEntryType == vtkFoamToken::SCALARLIST
            && entry.FirstValue().GetType() == vtkFoamToken::VECTORLIST
            && nRefComponents == nComponents))
          {
          for (int tupleI = 0; tupleI < nTuples; tupleI++)
            {
            float *tuple = data->GetPointer(nComponents * tupleI);
            for (int componentI = 0; componentI < nComponents; componentI++)
              {
              tuple[componentI] += refTuple[componentI];
              }
            }
          }
        else
          {
          vtkErrorMacro(<< "Wrong referenceLevel type");
          return NULL;
          }
        }

#if VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS
      // swap the components of symmTensor to match the component
      // names in paraview
      // OpenFOAM: XX, XY, XZ, YY, YZ, ZZ
      // ParaView: XX, YY, ZZ, XY, YZ, XZ
      // if VTK is sufficiently new then use SetComponentName() instead
      const int nComponents = data->GetNumberOfComponents();
      if (nComponents == 6)
        {
        for (int tupleI = 0; tupleI < nTuples; tupleI++)
          {
          float *tuple = data->GetPointer(nComponents * tupleI);
          const float symxy = tuple[1], symxz = tuple[2], symyy = tuple[3];
          const float symzz = tuple[5];
          tuple[1] = symyy;
          tuple[2] = symzz;
          tuple[3] = symxy;
          tuple[5] = symxz;
          }
        }
#endif
      }
    else if (entry.FirstValue().GetType() == vtkFoamToken::EMPTYLIST && nElements <= 0)
      {
      data = vtkFloatArray::New();
      // set the number of components as appropriate if the list is empty
      if (fieldType == "ScalarField" || fieldType == "SphericalTensorField")
        {
        data->SetNumberOfComponents(1);
        }
      else if (fieldType == "VectorField")
        {
        data->SetNumberOfComponents(3);
        }
      else if (fieldType == "SymmTensorField")
        {
        data->SetNumberOfComponents(6);
        }
      else if (fieldType == "TensorField")
        {
        data->SetNumberOfComponents(9);
        }
      }
    else
      {
      vtkErrorMacro(<< ioPtr->GetFileName().c_str() << " is not a valid "
          << ioPtr->GetClassName().c_str());
      return NULL;
      }
    }
  else // undefined
    {
    vtkErrorMacro(<< "Either uniform or nonuniform keyword is missing in entry "
        << entry.GetKeyword().c_str() << " in " << ioPtr->GetFileName().c_str());
    return NULL;
    }
  return data;
}

//-----------------------------------------------------------------------------
// convert OpenFOAM's dimension array representation to string
void vtkOFFReaderPrivate::ConstructDimensions(vtkStdString *dimString,
    vtkFoamDict *dictPtr)
{
  if (!this->Parent->GetAddDimensionsToArrayNames())
    {
    return;
    }
  vtkFoamEntry *dimEntry = dictPtr->Lookup("dimensions");
  if (dimEntry != NULL && dimEntry->FirstValue().GetType() == vtkFoamToken::LABELLIST)
    {
    vtkIntArray &dims = dimEntry->LabelList();
    if (dims.GetNumberOfTuples() == 7)
      {
      int dimSet[7];
      for (int dimI = 0; dimI < 7; dimI++)
        {
        dimSet[dimI] = dims.GetValue(dimI);
        }
      static const char *units[7] =
      { "kg", "m", "s", "K", "mol", "A", "cd" };
      vtksys_ios::ostringstream posDim, negDim;
      int posSpc = 0, negSpc = 0;
      if (dimSet[0] == 1 && dimSet[1] == -1 && dimSet[2] == -2)
        {
        posDim << "Pa";
        dimSet[0] = dimSet[1] = dimSet[2] = 0;
        posSpc = 1;
        }
      for (int dimI = 0; dimI < 7; dimI++)
        {
        const int dimDim = dimSet[dimI];
        if (dimDim > 0)
          {
          if (posSpc)
            {
            posDim << " ";
            }
          posDim << units[dimI];
          if (dimDim > 1)
            {
            posDim << dimDim;
            }
          posSpc++;
          }
        else if (dimDim < 0)
          {
          if (negSpc)
            {
            negDim << " ";
            }
          negDim << units[dimI];
          if (dimDim < -1)
            {
            negDim << -dimDim;
            }
          negSpc++;
          }
        }
      *dimString += " [" + posDim.str();
      if (negSpc > 0)
        {
        if (posSpc == 0)
          {
          *dimString += "1";
          }
        if (negSpc > 1)
          {
          *dimString += "/(" + negDim.str() + ")";
          }
        else
          {
          *dimString += "/" + negDim.str();
          }
        }
      else if (posSpc == 0)
        {
        *dimString += "-";
        }
      *dimString += "]";
      }
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::GetVolFieldAtTimeStep(
    vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
    const vtkStdString &varName)
{
  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkFoamDict dict;
  if (!this->ReadFieldFile(&io, &dict, varName,
      this->Parent->CellDataArraySelection))
    {
    return;
    }

  if (io.GetClassName().substr(0, 3) != "vol")
    {
    vtkErrorMacro(<< io.GetFileName().c_str() << " is not a volField");
    return;
    }

  const bool dimensionedInternalField
      = (io.GetClassName().find("::DimensionedInternalField")
      != vtkStdString::npos);
  vtkFoamEntry *iEntry
      = dict.Lookup(dimensionedInternalField ? "value" : "internalField");
  if (!iEntry)
    {
    vtkErrorMacro(<<"internalField not found in " << io.GetFileName().c_str());
    return;
    }

  if (iEntry->FirstValue().GetType() == vtkFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if (this->NumCells > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.GetFileName().c_str()
          << " is empty");
      }
    return;
    }

  vtkFoamEntry *rEntry = dict.Lookup("referenceLevel");
  vtkStdString fieldType = dimensionedInternalField
      ? io.GetClassName().substr(3, io.GetClassName().find_first_of(':') - 3)
      : io.GetClassName().substr(3, vtkStdString::npos);
  vtkFloatArray *iData =
      this->FillField(iEntry, rEntry, this->NumCells, &io, fieldType, false);
  if (iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  vtkFloatArray *ctpData = NULL;

  if (iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if (internalMesh != NULL)
      {
      if (this->Parent->GetDecomposePolyhedra())
        {
        // add values for decomposed cells
        this->ExtendArray<vtkFloatArray, float>(iData, this->NumCells
            + this->NumTotalAdditionalCells);
        const int nTuples = this->AdditionalCellIds->GetNumberOfTuples();
        int additionalCellI = this->NumCells;
        for (int tupleI = 0; tupleI < nTuples; tupleI++)
          {
          const int nCells = this->NumAdditionalCells->GetValue(tupleI);
          const vtkIdType cellId = this->AdditionalCellIds->GetValue(tupleI);
          for (int cellI = 0; cellI < nCells; cellI++)
            {
            iData->InsertTuple(additionalCellI++, cellId, iData);
            }
          }
        }

      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetCellData(), iData,
          io.GetObjectName() + dimString);

      if (this->Parent->GetCreateCellToPoint())
        {
        // Create cell-to-point interpolated data
        ctpData = vtkFloatArray::New();
        ctpData->SetNumberOfComponents(iData->GetNumberOfComponents());
        ctpData->SetNumberOfTuples(internalMesh->GetPoints()->GetNumberOfPoints());
        if (this->InternalPoints != NULL)
          {
          this->InterpolateCellToPoint(ctpData, iData, internalMesh,
              this->InternalPoints, this->InternalPoints->GetNumberOfTuples());
          }

        if (this->Parent->GetDecomposePolyhedra())
          {
          // assign cell values to additional points
          const int nPoints = this->AdditionalCellIds->GetNumberOfTuples();
          for (int pointI = 0; pointI < nPoints; pointI++)
            {
            ctpData->SetTuple(this->NumPoints + pointI,
                this->AdditionalCellIds->GetValue(pointI), iData);
            }
          }
        }
      }
    }
  else
    {
    // determine as there's no cells
    iData->Delete();
    return;
    }

  if (boundaryMesh == NULL || dimensionedInternalField)
    {
    iData->Delete();
    if (ctpData != NULL)
      {
      ctpData->Delete();
      }
    return;
    }

  vtkFloatArray *acData = NULL;

  if (this->Parent->GetCreateCellToPoint())
    {
    if (this->AllBoundaries == NULL)
      {
      vtkErrorMacro(<<"boundary mesh for cell to point filtering not found");
      iData->Delete();
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }
    acData = vtkFloatArray::New();
    acData->SetNumberOfComponents(iData->GetNumberOfComponents());
    acData->SetNumberOfTuples(this->AllBoundaries->GetNumberOfCells());
    }

  // set boundary values
  const vtkFoamEntry *bEntry = dict.Lookup("boundaryField");
  if (bEntry == NULL)
    {
    vtkErrorMacro(<< "boundaryField not found in object " << varName.c_str()
        << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
    iData->Delete();
    if (acData != NULL)
      {
      acData->Delete();
      }
    if (ctpData != NULL)
      {
      ctpData->Delete();
      }
    return;
    }

  for (int boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<int>(this->BoundaryDict.size()); boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const vtkStdString &boundaryNameI = beI.BoundaryName;

    const vtkFoamEntry *bEntryI = bEntry->Dictionary().Lookup(boundaryNameI);
    if (bEntryI == NULL)
      {
      vtkErrorMacro(<< "boundaryField " << boundaryNameI.c_str()
          << " not found in object " << varName.c_str() << " at time = "
          << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if (acData != NULL)
        {
        acData->Delete();
        }
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    if (bEntryI->FirstValue().GetType() != vtkFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<< "Type of boundaryField " << boundaryNameI.c_str()
          << " is not a subdictionary in object " << varName.c_str()
          << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if (acData != NULL)
        {
        acData->Delete();
        }
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    const int nFaces = beI.NFaces;

    vtkFloatArray* vData = NULL;
    if (beI.BoundaryType != vtkFoamBoundaryEntry::GEOMETRICAL)
      {
      vtkFoamEntry *vEntry = bEntryI->Dictionary().Lookup("value");
      if (vEntry != NULL) // the boundary has a value entry
        {
        vData = this->FillField(vEntry, rEntry, nFaces, &io, fieldType, false);
        if (vData == NULL)
          {
          iData->Delete();
          if (acData != NULL)
            {
            acData->Delete();
            }
          if (ctpData != NULL)
            {
            ctpData->Delete();
            }
          return;
          }
        }
      else
        {
        // uniformFixedValue B.C.
        const vtkFoamEntry *ufvEntry = bEntryI->Dictionary().Lookup("type");
        if (ufvEntry != NULL)
          {
          if (ufvEntry->ToWord() == "uniformFixedValue")
            {
            // the boundary is of uniformFixedValue type
            vtkFoamEntry *uvEntry = bEntryI->Dictionary().Lookup("uniformValue");
            if (uvEntry != NULL) // and has a uniformValue entry
              {
              vData = this->FillField(uvEntry, rEntry, nFaces, &io, fieldType, true);
              if (vData == NULL)
                {
                iData->Delete();
                if (acData != NULL)
                  {
                  acData->Delete();
                  }
                if (ctpData != NULL)
                  {
                  ctpData->Delete();
                  }
                return;
                }
              }
            }
          }
        }
      }
    const int boundaryStartFace = beI.StartFace
        - this->BoundaryDict[0].StartFace;

    if (vData == NULL) // doesn't have a value nor uniformValue entry
      {
      // use patch-internal values as boundary values
      vData = vtkFloatArray::New();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nFaces);
      for (int j = 0; j < nFaces; j++)
        {
        const int cellId = this->FaceOwner->GetValue(boundaryStartFace + j);
        vData->SetTuple(j, cellId, iData);
        }
      }

    if (this->Parent->GetCreateCellToPoint())
      {
      const int startFace = beI.AllBoundariesStartFace;
      // if reading a processor sub-case of a decomposed case as is,
      // use the patch values of the processor patch as is
      if (beI.BoundaryType == vtkFoamBoundaryEntry::PHYSICAL
          || (this->ProcessorName == "" && beI.BoundaryType
              == vtkFoamBoundaryEntry::PROCESSOR))
        {
        // set the same value to AllBoundaries
        for (int faceI = 0; faceI < nFaces; faceI++)
          {
          acData->SetTuple(faceI + startFace, faceI, vData);
          }
        }
      // implies && this->ProcessorName != ""
      else if (beI.BoundaryType == vtkFoamBoundaryEntry::PROCESSOR)
        {
        // average patch internal value and patch value assuming the
        // patch value to be the patchInternalField of the neighbor
        // decomposed mesh. Using double precision to avoid degrade in
        // accuracy.
        const int nComponents = vData->GetNumberOfComponents();
        for (int faceI = 0; faceI < nFaces; faceI++)
          {
          const float *vTuple = vData->GetPointer(nComponents * faceI);
          const float *iTuple = iData->GetPointer(nComponents
              * this->FaceOwner->GetValue(boundaryStartFace + faceI));
          float *acTuple =
              acData->GetPointer(nComponents * (startFace + faceI));
          for (int componentI = 0; componentI < nComponents; componentI++)
            {
            acTuple[componentI] = (static_cast<double>(vTuple[componentI])
                + static_cast<double>(iTuple[componentI])) * 0.5;
            }
          }
        }
      }

    if (beI.IsActive)
      {
      vtkPolyData *bm =
          vtkPolyData::SafeDownCast(boundaryMesh->GetBlock(activeBoundaryI));
      this->AddArrayToFieldData(bm->GetCellData(), vData, io.GetObjectName()
          + dimString);

      if (this->Parent->GetCreateCellToPoint())
        {
        // construct cell-to-point interpolated boundary values. This
        // is done independently from allBoundary interpolation so
        // that the interpolated values are not affected by
        // neighboring patches especially at patch edges and for
        // baffle patches
        vtkFloatArray *pData = vtkFloatArray::New();
        pData->SetNumberOfComponents(vData->GetNumberOfComponents());
        const int nPoints = bm->GetPoints()->GetNumberOfPoints();
        pData->SetNumberOfTuples(nPoints);
        this->InterpolateCellToPoint(pData, vData, bm, NULL, nPoints);
        this->AddArrayToFieldData(bm->GetPointData(), pData, io.GetObjectName()
            + dimString);
        pData->Delete();
        }

      activeBoundaryI++;
      }
    vData->Delete();
    }
  iData->Delete();

  if (this->Parent->GetCreateCellToPoint())
    {
    // Create cell-to-point interpolated data for all boundaries and
    // override internal values
    vtkFloatArray *bpData = vtkFloatArray::New();
    bpData->SetNumberOfComponents(acData->GetNumberOfComponents());
    const int nPoints = this->AllBoundariesPointMap->GetNumberOfTuples();
    bpData->SetNumberOfTuples(nPoints);
    this->InterpolateCellToPoint(bpData, acData, this->AllBoundaries, NULL,
        nPoints);
    acData->Delete();

    if (ctpData != NULL)
      {
      // set cell-to-pint data for internal mesh
      for (int pointI = 0; pointI < nPoints; pointI++)
        {
        ctpData->SetTuple(this->AllBoundariesPointMap->GetValue(pointI),
            pointI, bpData);
        }
      this->AddArrayToFieldData(internalMesh->GetPointData(), ctpData,
          io.GetObjectName() + dimString);
      ctpData->Delete();
      }

    bpData->Delete();
    }
}

//-----------------------------------------------------------------------------
// read surface field at a timestep
void vtkOFFReaderPrivate::GetSurfaceFieldAtTimeStep(
    vtkPolyData *surfaceMesh, vtkMultiBlockDataSet *boundaryMesh,
    const vtkStdString &varName)
{
  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkFoamDict dict;
  if (!this->ReadFieldFile(&io, &dict, varName,
      this->Parent->SurfaceDataArraySelection))
    {
    return;
    }

  if (io.GetClassName().substr(0, 7) != "surface")
    {
    vtkErrorMacro(<< io.GetFileName().c_str() << " is not a surfaceField");
    return;
    }

  vtkFoamEntry *iEntry = dict.Lookup("internalField");
  if (iEntry == NULL)
    {
    vtkErrorMacro (<< "internalField not found in "
        << io.GetFileName().c_str());
    return;
    }

  vtkFoamEntry *rEntry = dict.Lookup("referenceLevel");
  vtkStdString fieldType = (io.GetClassName().substr(7, vtkStdString::npos));

  vtkFloatArray *iData = NULL;
  if (surfaceMesh != NULL)
    {
    // Check if internal faces have been defined
    if (this->NumInternalFaces <= 0)
      {
      vtkErrorMacro(<< "Internal faces is either invalid or not defined."
          << " Internal Faces: " << this->NumInternalFaces);
      return;
      }

    iData = this->FillField(iEntry, rEntry, this->NumInternalFaces, &io,
        fieldType, false);
    if (iData == NULL)
      {
      return;
      }

    if (this->ProcessorFaces != NULL)
      {
      // the extended array will be filled later
      this->ExtendArray<vtkFloatArray, float>(iData, this->NumInternalFaces
          + this->ProcessorFaces->GetBody()->GetNumberOfTuples());
      }
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  if (boundaryMesh != NULL || this->ProcessorFaces != NULL)
    {
    // set boundary values
    const vtkFoamEntry *bEntry = dict.Lookup("boundaryField");
    if (bEntry == NULL)
      {
      vtkErrorMacro( << "boundaryField not found in object " << varName.c_str()
          << " at time = "
          << this->TimeNames->GetValue(this->TimeStep).c_str());
      if (iData != NULL)
        {
        iData->Delete();
        }
      return;
      }

    for (int boundaryI = 0, activeBoundaryI = 0, procBoundaryI = 0;
        boundaryI < static_cast<int>(this->BoundaryDict.size()); boundaryI++)
      {
      const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];

      const bool doProcFaces = (this->ProcessorFaces != NULL && procBoundaryI
          < static_cast<int>(this->BoundaryDict.ProcBoundaries.size())
          && this->BoundaryDict.ProcBoundaries[procBoundaryI] == boundaryI);
      if (!beI.IsActive && !doProcFaces)
        {
        continue;
        }

      const vtkStdString &boundaryNameI = beI.BoundaryName;
      const vtkFoamEntry *bEntryI = bEntry->Dictionary().Lookup(boundaryNameI);

      if (bEntryI == NULL)
        {
        vtkErrorMacro(<< "boundaryField " << boundaryNameI.c_str()
            << " not found in object " << varName.c_str() << " at time = "
            << this->TimeNames->GetValue(this->TimeStep).c_str());

        if (iData != NULL)
          {
          iData->Delete();
          }
        return;
        }

      if (bEntryI->FirstValue().GetType() != vtkFoamToken::DICTIONARY)
        {
        vtkErrorMacro( << "Type of boundaryField " << boundaryNameI.c_str()
            << " is not a subdictionary in object " << varName.c_str()
            << " at time = "
            << this->TimeNames->GetValue(this->TimeStep).c_str());

        if (iData != NULL)
          {
          iData->Delete();
          }
        return;
        }

      const int nFaces = beI.NFaces;
      vtkFloatArray* vData = NULL;

      if (beI.BoundaryType != vtkFoamBoundaryEntry::GEOMETRICAL)
        {
        vtkFoamEntry *vEntry = bEntryI->Dictionary().Lookup("value");

        if (vEntry != NULL)
          {
          // The boundary has a value entry
          vData = this->FillField(vEntry, rEntry, nFaces, &io, fieldType,
              false);
          if (vData == NULL)
            {
            if (iData != NULL)
              {
              iData->Delete();
              }
            return;
            }
          }
        else
          {
          // uniformFixedValue B.C.
          const vtkFoamEntry *ufvEntry = bEntryI->Dictionary().Lookup("type");
          if (ufvEntry != NULL)
            {
            if (ufvEntry->ToWord() == "uniformFixedValue")
              {
              // the boundary is of uniformFixedValue type
              vtkFoamEntry *uvEntry
                = bEntryI->Dictionary().Lookup("uniformValue");
              if (uvEntry != NULL)
                {
                // and has a uniformValue entry
                vData = this->FillField(uvEntry, rEntry, nFaces, &io, fieldType,
                    true);
                if (vData == NULL)
                  {
                  if (iData != NULL)
                    {
                    iData->Delete();
                    }
                  return;
                  }
                }
              }
            }
          }
        }

      if (vData == NULL)
        {
        // Doesn't have a value nor uniformValue entry.
        // Cannot do the equivalent of zero-gradient in cell-fields,
        // so this will be artificially set to zero.
        vData = vtkFloatArray::New();
        if (fieldType == "ScalarField" || fieldType == "SphericalTensorField")
          {
          vData->SetNumberOfComponents(1);
          }
        else if (fieldType == "VectorField")
          {
          vData->SetNumberOfComponents(3);
          }
        else if (fieldType == "SymmTensorField")
          {
          vData->SetNumberOfComponents(6);
          }
        else if (fieldType == "TensorField")
          {
          vData->SetNumberOfComponents(9);
          }
        else
          {
          vtkErrorMacro(<< io.GetFileName().c_str() << " is not a valid field");
          if (iData != NULL)
            {
            iData->Delete();
            }
          return;
          }
        vData->SetNumberOfTuples(nFaces);

        for (int j = 0; j < vData->GetNumberOfComponents(); j++)
          {
          vData->FillComponent(j, 0.0);
          }
        }

      if (beI.IsActive)
        {
        vtkPolyData *bm = vtkPolyData::SafeDownCast(
            boundaryMesh->GetBlock(activeBoundaryI));
        this->AddArrayToFieldData(bm->GetCellData(), vData,
            io.GetObjectName() + dimString);
        activeBoundaryI++;
        }

      if (doProcFaces)
        {
        const int startFace = this->NumInternalFaces
            + this->ProcessorFaces->GetIndices()->GetValue(procBoundaryI);
        const int nProcFaces = this->ProcessorFaces->GetSize(procBoundaryI);
        const int *faceIDs = this->ProcessorFaces->operator[](procBoundaryI);
        for (int faceI = 0; faceI < nProcFaces; faceI++)
          {
          iData->SetTuple(startFace + faceI, faceIDs[faceI], vData);
          }
        procBoundaryI++;
        }
      vData->Delete();
      }
    }

  // Skip surface mesh, if it's not requested
  if (surfaceMesh != NULL)
    {
    if (iData->GetSize() > 0)
      {
      // set data to internal mesh
      this->AddArrayToFieldData(surfaceMesh->GetCellData(), iData,
          io.GetObjectName() + dimString);
      }
    iData->Delete();
    }
}

//-----------------------------------------------------------------------------
// read point field at a timestep
void vtkOFFReaderPrivate::GetPointFieldAtTimeStep(
    vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
    const vtkStdString &varName)
{
  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkFoamDict dict;
  if (!this->ReadFieldFile(&io, &dict, varName,
      this->Parent->PointDataArraySelection))
    {
    return;
    }

  if (io.GetClassName().substr(0, 5) != "point")
    {
    vtkErrorMacro(<< io.GetFileName().c_str() << " is not a pointField");
    return;
    }

  vtkFoamEntry *iEntry = dict.Lookup("internalField");
  if (iEntry == NULL)
    {
    vtkErrorMacro(<<"internalField not found in " << io.GetFileName().c_str());
    return;
    }

  if (iEntry->FirstValue().GetType() == vtkFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if (this->NumPoints > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.GetFileName().c_str()
          << " is empty");
      }
    return;
    }

  vtkFoamEntry *rEntry = dict.Lookup("referenceLevel");
  vtkStdString fieldType = io.GetClassName().substr(5, vtkStdString::npos);
  vtkFloatArray *iData = this->FillField(iEntry, rEntry, this->NumPoints, &io,
      fieldType, false);
  if (iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  // AdditionalCellPoints is NULL if creation of InternalMesh had been skipped
  if (this->AdditionalCellPoints != NULL)
    {
    // point-to-cell interpolation to additional cell centroidal points
    // for decomposed cells
    const int nAdditionalPoints = static_cast<int>(this->AdditionalCellPoints->size());
    const int nComponents = iData->GetNumberOfComponents();
    this->ExtendArray<vtkFloatArray, float>(iData, this->NumPoints
        + nAdditionalPoints);
    for (int i = 0; i < nAdditionalPoints; i++)
      {
      vtkIdList *acp = this->AdditionalCellPoints->operator[](i);
      const vtkIdType nPoints = acp->GetNumberOfIds();
      double interpolatedValue[9];
      for (int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] = 0.0;
        }
      for (vtkIdType j = 0; j < nPoints; j++)
        {
        const float *tuple = iData->GetPointer(nComponents * acp->GetId(j));
        for (int k = 0; k < nComponents; k++)
          {
          interpolatedValue[k] += tuple[k];
          }
        }
      const double weight = 1.0 / static_cast<double>(nPoints);
      for (int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] *= weight;
        }
      // will automatically be converted to float
      iData->InsertTuple(this->NumPoints + i, interpolatedValue);
      }
    }

  if (iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if (internalMesh != NULL)
      {
      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetPointData(), iData,
          io.GetObjectName() + dimString);
      }
    }
  else
    {
    // determine as there's no points
    iData->Delete();
    return;
    }

  if (boundaryMesh == NULL)
    {
    iData->Delete();
    return;
    }

  // use patch-internal values as boundary values
  for (vtkIdType boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<vtkIdType>(this->BoundaryDict.size()); boundaryI++)
    {
    if (this->BoundaryDict[boundaryI].IsActive)
      {
      vtkFloatArray *vData = vtkFloatArray::New();
      vtkIntArray& bpMap = *this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nPoints = bpMap.GetNumberOfTuples();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nPoints);
      for (int j = 0; j < nPoints; j++)
        {
        vData->SetTuple(j, bpMap.GetValue(j), iData);
        }
      this->AddArrayToFieldData(vtkPolyData::SafeDownCast(
          boundaryMesh->GetBlock(activeBoundaryI))->GetPointData(), vData, io.GetObjectName()
          + dimString);
      vData->Delete();
      activeBoundaryI++;
      }
    }
  iData->Delete();
}

//-----------------------------------------------------------------------------
vtkMultiBlockDataSet* vtkOFFReaderPrivate::MakeLagrangianMesh()
{
  vtkMultiBlockDataSet *lagrangianMesh = vtkMultiBlockDataSet::New();

  for (int cloudI = 0; cloudI
      < this->Parent->LagrangianPaths->GetNumberOfTuples(); cloudI++)
    {
    const vtkStdString& pathI = this->Parent->LagrangianPaths->GetValue(cloudI);

    // still can't distinguish on patch selection panel, but can
    // distinguish the "lagrangian" reserved path component and a mesh
    // region with the same name
    vtkStdString subCloudName;
    if (pathI[0] == '/')
      {
      subCloudName = pathI.substr(1, vtkStdString::npos);
      }
    else
      {
      subCloudName = pathI;
      }
    if (this->RegionName != pathI.substr(0, pathI.find('/'))
        || !this->Parent->GetPatchArrayStatus(subCloudName.c_str()))
      {
      continue;
      }

    const vtkStdString cloudPath(this->CurrentTimePath() + "/" + subCloudName
        + "/");
    const vtkStdString positionsPath(cloudPath + "positions");

    // create an empty mesh to keep node/leaf structure of the
    // multi-block consistent even if mesh doesn't exist
    vtkPolyData *meshI = vtkPolyData::New();
    const int blockI = lagrangianMesh->GetNumberOfBlocks();
    lagrangianMesh->SetBlock(blockI, meshI);
    // extract the cloud name
    this->SetBlockName(lagrangianMesh, blockI, pathI.substr(pathI.rfind('/') + 1).c_str());

    vtkFoamIOobject io(this->CasePath,
        this->Parent->GetIsSinglePrecisionBinary() != 0);
    if (!(io.Open(positionsPath) || io.Open(positionsPath + ".gz")))
      {
      meshI->Delete();
      continue;
      }

    // tell the IO object if the file is in OF 1.3 binary
    // lagrangian/positions format
    io.SetIs13Positions(this->Parent->GetPositionsIsIn13Format() != 0);

    vtkFoamEntryValue dict(NULL);
    try
      {
      dict.ReadNonuniformList<vtkFoamToken::VECTORLIST,
      vtkFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, true> >(
          io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      meshI->Delete();
      continue;
      }
    io.Close();

    vtkFloatArray *pointArray = static_cast<vtkFloatArray *>(dict.Ptr());
    const int nParticles = pointArray->GetNumberOfTuples();

    // instantiate the points class
    vtkPoints *points = vtkPoints::New();
    points->SetData(pointArray);
    pointArray->Delete();

    // create lagrangian mesh
    meshI->Allocate(nParticles);
    for (vtkIdType i = 0; i < nParticles; i++)
      {
      meshI->InsertNextCell(VTK_VERTEX, 1, &i);
      }
    meshI->SetPoints(points);
    points->Delete();

    // read lagrangian fields
    for (int fieldI = 0; fieldI
        < this->LagrangianFieldFiles->GetNumberOfValues(); fieldI++)
      {
      const vtkStdString varPath(cloudPath
          + this->LagrangianFieldFiles->GetValue(fieldI));

      vtkFoamIOobject io2(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      if (!io2.Open(varPath))
        {
        // if the field file doesn't exist we simply return without
        // issuing an error as a simple way of supporting multi-region
        // lagrangians
        continue;
        }

      // if the variable is disabled on selection panel then skip it
      const vtkStdString selectionName(io2.GetObjectName());
      if (this->Parent->LagrangianDataArraySelection->ArrayExists(selectionName.c_str())
          && !this->Parent->GetLagrangianArrayStatus(selectionName.c_str()))
        {
        continue;
        }

      // read the field file into dictionary
      vtkFoamEntryValue dict2(NULL);
      if (!dict2.ReadField(io2))
        {
        vtkErrorMacro(<<"Error reading line " << io2.GetLineNumber()
            << " of " << io2.GetFileName().c_str() << ": "
            << io2.GetError().c_str());
        continue;
        }

      // set lagrangian values
      if (dict2.GetType() != vtkFoamToken::SCALARLIST && dict2.GetType()
          != vtkFoamToken::VECTORLIST && dict2.GetType()
          != vtkFoamToken::LABELLIST)
        {
        vtkErrorMacro(<< io2.GetFileName().c_str()
            << ": Unsupported lagrangian field type "
            << io2.GetClassName().c_str());
        continue;
        }

      vtkDataArray* lData = static_cast<vtkDataArray *>(dict2.Ptr());

      // GetNumberOfTuples() works for both scalar and vector
      const int nParticles2 = lData->GetNumberOfTuples();
      if (nParticles2 != meshI->GetNumberOfCells())
        {
        vtkErrorMacro(<< io2.GetFileName().c_str()
            <<": Sizes of lagrangian mesh and field don't match: mesh = "
            << meshI->GetNumberOfCells() << ", field = " << nParticles2);
        lData->Delete();
        continue;
        }

      this->AddArrayToFieldData(meshI->GetCellData(), lData, selectionName);
      if (this->Parent->GetCreateCellToPoint())
        {
        this->AddArrayToFieldData(meshI->GetPointData(), lData, selectionName);
        }
      lData->Delete();
      }
    meshI->Delete();
    }
  return lagrangianMesh;
}

//-----------------------------------------------------------------------------
// returns a dictionary of block names for a specified domain
vtkFoamDict* vtkOFFReaderPrivate::GatherBlocks(const char* typeIn,
    const int timeStep)
{
  if (this->PolyMeshFacesDir->GetValue(timeStep) == "")
    {
    return NULL;
    }

  vtkStdString type(typeIn);
  vtkStdString blockPath =
      this->TimeRegionMeshPath(this->PolyMeshFacesDir, timeStep) + type;

  vtkFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(blockPath) || io.Open(blockPath + ".gz")))
    {
    return NULL;
    }

  vtkFoamDict* dictPtr = new vtkFoamDict;
  vtkFoamDict& dict = *dictPtr;
  if (!dict.Read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
    delete dictPtr;
    return NULL;
    }
  if (dict.GetType() != vtkFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"The file type of " << io.GetFileName().c_str()
        << " is not a dictionary");
    delete dictPtr;
    return NULL;
    }
  return dictPtr;
}

//-----------------------------------------------------------------------------
// returns a requested point zone mesh
bool vtkOFFReaderPrivate::GetPointZoneMesh(
    vtkMultiBlockDataSet *pointZoneMesh, vtkPoints *points)
{
  vtkFoamDict *pointZoneDictPtr
      = this->GatherBlocks("pointZones", this->TimeStep);

  if (pointZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &pointZoneDict = *pointZoneDictPtr;
  int nPointZones = static_cast<int>(pointZoneDict.size());

  for (int i = 0; i < nPointZones; i++)
    {
    // look up point labels
    vtkFoamDict &dict = pointZoneDict[i]->Dictionary();
    vtkFoamEntry *pointLabelsEntry = dict.Lookup("pointLabels");
    if (pointLabelsEntry == NULL)
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not found in pointZones");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (pointLabelsEntry->FirstValue().GetType() == vtkFoamToken::EMPTYLIST)
      {
      vtkPolyData *pzm = vtkPolyData::New();
      pointZoneMesh->SetBlock(i, pzm);
      pzm->Delete();
      // set name
      this->SetBlockName(pointZoneMesh, i, pointZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (pointLabelsEntry->FirstValue().GetType() != vtkFoamToken::LABELLIST)
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not of type labelList: type = "
          << pointLabelsEntry->FirstValue().GetType());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    vtkIntArray &labels = pointLabelsEntry->LabelList();

    int nPoints = labels.GetNumberOfTuples();
    if (nPoints > this->NumPoints)
      {
      vtkErrorMacro(<<"The length of pointLabels " << nPoints
          << " for pointZone " << pointZoneDict[i]->GetKeyword().c_str()
          << " exceeds the number of points " << this->NumPoints);
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      delete pointZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkPolyData *pzm = vtkPolyData::New();

    // set pointZone size
    pzm->Allocate(nPoints);

    // insert points
    for (int j = 0; j < nPoints; j++)
      {
      vtkIdType pointLabel = labels.GetValue(j); // must be vtkIdType
      if (pointLabel >= this->NumPoints)
        {
        vtkWarningMacro(<<"pointLabels id " << pointLabel
            << " exceeds the number of points " << this->NumPoints);
        pzm->InsertNextCell(VTK_EMPTY_CELL, 0, &pointLabel);
        continue;
        }
      pzm->InsertNextCell(VTK_VERTEX, 1, &pointLabel);
      }
    pzm->SetPoints(points);

    pointZoneMesh->SetBlock(i, pzm);
    pzm->Delete();
    // set name
    this->SetBlockName(pointZoneMesh, i, pointZoneDict[i]->GetKeyword().c_str());
    }

  delete pointZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested face zone mesh
bool vtkOFFReaderPrivate::GetFaceZoneMesh(
    vtkMultiBlockDataSet *faceZoneMesh,
    const vtkFoamIntVectorVector *facesPoints, vtkPoints *points)
{
  vtkFoamDict *faceZoneDictPtr
      = this->GatherBlocks("faceZones", this->TimeStep);

  if (faceZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &faceZoneDict = *faceZoneDictPtr;
  int nFaceZones = static_cast<int>(faceZoneDict.size());

  for (int i = 0; i < nFaceZones; i++)
    {
    // look up face labels
    vtkFoamDict &dict = faceZoneDict[i]->Dictionary();
    vtkFoamEntry *faceLabelsEntry = dict.Lookup("faceLabels");
    if (faceLabelsEntry == NULL)
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not found in faceZones");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (faceLabelsEntry->FirstValue().GetType() == vtkFoamToken::EMPTYLIST)
      {
      vtkPolyData *fzm = vtkPolyData::New();
      faceZoneMesh->SetBlock(i, fzm);
      fzm->Delete();
      // set name
      this->SetBlockName(faceZoneMesh, i, faceZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (faceLabelsEntry->FirstValue().GetType() != vtkFoamToken::LABELLIST)
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not of type labelList");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    vtkIntArray &labels = faceLabelsEntry->LabelList();

    int nFaces = labels.GetNumberOfTuples();
    if (nFaces > this->FaceOwner->GetNumberOfTuples())
      {
      vtkErrorMacro(<<"The length of faceLabels " << nFaces
          << " for faceZone " << faceZoneDict[i]->GetKeyword().c_str()
          << " exceeds the number of faces "
          << this->FaceOwner->GetNumberOfTuples());
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      delete faceZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkPolyData *fzm = vtkPolyData::New();

    // set faceZone size
    fzm->Allocate(nFaces);

    // aloocate array for converting int vector to vtkIdType vector:
    // workaround for 64bit machines
    int maxNFacePoints = 0;
    for (int j = 0; j < nFaces; j++)
      {
      const int nFacePoints = facesPoints->GetSize(labels.GetValue(j));
      if (nFacePoints > maxNFacePoints)
        {
        maxNFacePoints = nFacePoints;
        }
      }
    vtkIdList *facePointsVtkId = vtkIdList::New();
    facePointsVtkId->SetNumberOfIds(maxNFacePoints);

    // insert faces
    if (!this->InsertFacesToGrid(fzm, facesPoints, 0, nFaces, NULL,
        facePointsVtkId, labels.GetPointer(0), false))
      {
      facePointsVtkId->Delete();
      fzm->Delete();
      return false;
      }
    facePointsVtkId->Delete();
    fzm->SetPoints(points);
    faceZoneMesh->SetBlock(i, fzm);
    fzm->Delete();
    // set name
    this->SetBlockName(faceZoneMesh, i, faceZoneDict[i]->GetKeyword().c_str());
    }

  delete faceZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested cell zone mesh
bool vtkOFFReaderPrivate::GetCellZoneMesh(
    vtkMultiBlockDataSet *cellZoneMesh,
    const vtkFoamIntVectorVector *cellsFaces,
    const vtkFoamIntVectorVector *facesPoints, vtkPoints *points)
{
  vtkFoamDict *cellZoneDictPtr
      = this->GatherBlocks("cellZones", this->TimeStep);

  if (cellZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &cellZoneDict = *cellZoneDictPtr;
  int nCellZones = static_cast<int>(cellZoneDict.size());

  for (int i = 0; i < nCellZones; i++)
    {
    // look up cell labels
    vtkFoamDict &dict = cellZoneDict[i]->Dictionary();
    vtkFoamEntry *cellLabelsEntry = dict.Lookup("cellLabels");
    if (cellLabelsEntry == NULL)
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not found in cellZones");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (cellLabelsEntry->FirstValue().GetType() == vtkFoamToken::EMPTYLIST)
      {
      vtkUnstructuredGrid *czm = vtkUnstructuredGrid::New();
      cellZoneMesh->SetBlock(i, czm);
      // set name
      this->SetBlockName(cellZoneMesh, i, cellZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (cellLabelsEntry->FirstValue().GetType() != vtkFoamToken::LABELLIST)
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not of type labelList");
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      return false;
      }

    vtkIntArray &labels = cellLabelsEntry->LabelList();

    int nCells = labels.GetNumberOfTuples();
    if (nCells > this->NumCells)
      {
      vtkErrorMacro(<<"The length of cellLabels " << nCells
          << " for cellZone " << cellZoneDict[i]->GetKeyword().c_str()
          << " exceeds the number of cells " << this->NumCells);
      this->Parent->SetErrorCode(vtkErrorCode::FileFormatError);
      delete cellZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointers if we return by error
    vtkUnstructuredGrid *czm = vtkUnstructuredGrid::New();

    // set cellZone size
    czm->Allocate(nCells);

    // insert cells
    this->InsertCellsToGrid(czm, cellsFaces, facesPoints, NULL, NULL, &labels);

    // set cell zone points
    czm->SetPoints(points);

    cellZoneMesh->SetBlock(i, czm);
    czm->Delete();

    // set name
    this->SetBlockName(cellZoneMesh, i, cellZoneDict[i]->GetKeyword().c_str());
    }

  delete cellZoneDictPtr;
  return true;
}

//-----------------------------------------------------------------------------
void vtkOFFReaderPrivate::AddArrayToFieldData(
    vtkDataSetAttributes *fieldData, vtkDataArray *array,
    const vtkStdString &arrayName)
{
  // exclude dimensional unit string if any
  const vtkStdString arrayNameString(arrayName.substr(0, arrayName.find(' ')));
  array->SetName(arrayName.c_str());

#if !VTK_FOAMFILE_REORDER_SYMMTENSOR_COMPONENTS
  switch (array->GetNumberOfComponents())
    {
    case 3:
      array->SetComponentName(0, "x");
      array->SetComponentName(1, "y");
      array->SetComponentName(2, "z");
      break;
    case 6:
      array->SetComponentName(0, "xx");
      array->SetComponentName(1, "xy");
      array->SetComponentName(2, "xz");
      array->SetComponentName(3, "yy");
      array->SetComponentName(4, "yz");
      array->SetComponentName(5, "zz");
      break;
    case 9:
      array->SetComponentName(0, "xx");
      array->SetComponentName(1, "xy");
      array->SetComponentName(2, "xz");
      array->SetComponentName(3, "yx");
      array->SetComponentName(4, "yy");
      array->SetComponentName(5, "yz");
      array->SetComponentName(6, "zx");
      array->SetComponentName(7, "zy");
      array->SetComponentName(8, "zz");
      break;
    }
#endif

  if (array->GetNumberOfComponents() == 1 && arrayNameString == "p")
    {
    fieldData->SetScalars(array);
    }
  else if (array->GetNumberOfComponents() == 3 && arrayNameString == "U")
    {
    fieldData->SetVectors(array);
    }
  else
    {
    fieldData->AddArray(array);
    }
}

//-----------------------------------------------------------------------------
// return 0 if there's any error, 1 if success
int vtkOFFReaderPrivate::RequestData(vtkMultiBlockDataSet *output,
bool recreateInternalMesh, bool recreateBoundaryMesh, bool updateVariables,
bool recreateLagrangianMesh)
{
  // determine if we need to reconstruct meshes
  recreateInternalMesh |= this->TimeStepOld == -1
      // the following three quite likely indicate reading mesh failed on
      // the previous RequestData() call, hence trying again
      || (this->InternalMeshSelectionStatus && this->InternalMesh == NULL)
      || (this->SurfaceMeshSelectionStatus && this->SurfaceMesh == NULL)
      || this->BoundaryMesh == NULL
      || this->InternalMeshSelectionStatus
          != this->InternalMeshSelectionStatusOld
      || this->SurfaceMeshSelectionStatus
          != this->SurfaceMeshSelectionStatusOld
      || this->PolyMeshFacesDir->GetValue(this->TimeStep)
          != this->PolyMeshFacesDir->GetValue(this->TimeStepOld)
      || this->FaceOwner == NULL;
  recreateBoundaryMesh |= recreateInternalMesh;
  updateVariables |= recreateBoundaryMesh || this->TimeStep
      != this->TimeStepOld;
  recreateLagrangianMesh |= this->TimeStepOld == -1
      || (this->Parent->LagrangianPaths->GetNumberOfTuples()
          && this->TimeStep != this->TimeStepOld);
  const bool pointsMoved = this->TimeStepOld == -1
      || this->PolyMeshPointsDir->GetValue(this->TimeStep)
          != this->PolyMeshPointsDir->GetValue(this->TimeStepOld);
  const bool moveInternalPoints = !recreateInternalMesh && pointsMoved;
  const bool moveBoundaryPoints = !recreateBoundaryMesh && pointsMoved;

  const bool createEulerians
      = this->PolyMeshFacesDir->GetValue(this->TimeStep) != "";

  this->Parent->MeshChanged
      |= recreateBoundaryMesh || recreateLagrangianMesh || pointsMoved;

  if (recreateInternalMesh)
    {
    this->ClearInternalMeshes();
    }
  if (recreateBoundaryMesh)
    {
    this->ClearBoundaryMeshes();
    }
  if (recreateLagrangianMesh)
    {
    this->ClearLagrangianMeshes();
    }

  vtkFoamIntVectorVector *facePoints = NULL;
  vtkStdString meshDir;
  if (createEulerians && (recreateInternalMesh || recreateBoundaryMesh))
    {
    // create paths to polyMesh files
    meshDir = this->CurrentTimeRegionMeshPath(this->PolyMeshFacesDir);

    // create the faces vector
    facePoints = this->ReadFacesFile(meshDir);
    if (facePoints == NULL)
      {
      return 0;
      }
    this->Parent->UpdateProgress(0.2);
    }

  vtkFoamIntVectorVector *cellFaces = NULL;
  if (createEulerians && recreateInternalMesh)
    {
    // read owner/neighbor and create the FaceOwner and cellFaces vectors
    cellFaces = this->ReadOwnerNeighborFiles(meshDir, facePoints);
    if (cellFaces == NULL)
      {
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.3);
    }

  vtkFloatArray *pointArray = NULL;
  if (createEulerians && (recreateInternalMesh || (recreateBoundaryMesh
      && !recreateInternalMesh && this->InternalMesh == NULL)
      || moveInternalPoints || moveBoundaryPoints))
    {
    // get the points
    pointArray = this->ReadPointsFile();
    if ((pointArray == NULL && recreateInternalMesh) || (facePoints != NULL
        && !this->CheckFacePoints(facePoints)))
      {
      delete cellFaces;
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.4);
    }

  // make internal mesh
  // Create Internal Mesh only if required for display
  if (createEulerians && recreateInternalMesh)
    {
    if (this->InternalMeshSelectionStatus)
      {
      this->InternalMesh = this->MakeInternalMesh(cellFaces, facePoints,
          pointArray);
      }

    // Create the surface mesh, if necessary
    if (this->SurfaceMeshSelectionStatus)
      {
      vtkFloatArray *points;

      if (this->InternalMesh != NULL)
        {
        points = vtkFloatArray::SafeDownCast(
            this->InternalMesh->GetPoints()->GetData());
        }
      else
        {
        points = pointArray;
        }

      this->SurfaceMesh = this->MakeSurfaceMesh(meshDir, facePoints, points);

      if (this->SurfaceMesh == NULL)
        {
        delete cellFaces;
        delete facePoints;

        if (pointArray != NULL)
          {
          pointArray->Delete();
          }

        return 0;
        }
      }

    // read and construct zones
    if (this->Parent->GetReadZones())
      {
      vtkPoints *points;
      if (this->InternalMesh != NULL)
        {
        points = this->InternalMesh->GetPoints();
        }
      else
        {
        points = vtkPoints::New();
        points->SetData(pointArray);
        }

      this->PointZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetPointZoneMesh(this->PointZoneMesh, points))
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->PointZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        }

      this->FaceZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetFaceZoneMesh(this->FaceZoneMesh, facePoints, points))
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        if (this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->FaceZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        }

      this->CellZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetCellZoneMesh(this->CellZoneMesh, cellFaces, facePoints,
          points))
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        if (this->FaceZoneMesh != NULL)
          {
          this->FaceZoneMesh->Delete();
          this->FaceZoneMesh = NULL;
          }
        if (this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->CellZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        }
      if (this->InternalMesh == NULL)
        {
        points->Delete();
        }
      }
    delete cellFaces;
    this->TruncateFaceOwner();
    }

  if (createEulerians && recreateBoundaryMesh)
    {
    vtkFloatArray *boundaryPointArray;
    if (pointArray != NULL)
      {
      boundaryPointArray = pointArray;
      }
    else
      {
      boundaryPointArray
          = static_cast<vtkFloatArray *>(this->InternalMesh->GetPoints()->GetData());
      }
    // create boundary mesh
    this->BoundaryMesh = this->MakeBoundaryMesh(facePoints, boundaryPointArray);
    if (this->BoundaryMesh == NULL)
      {
      delete facePoints;
      if (pointArray != NULL)
        {
        pointArray->Delete();
        }
      return 0;
      }
#if 0
    this->CalculateReciprocalDelta(facePoints);
#endif
    }

  delete facePoints;

  // if only point coordinates change refresh point vector
  if (createEulerians && moveInternalPoints)
    {
    // refresh the points in each mesh
    vtkPoints *points;
    // Check if Internal Mesh exists first....
    if (this->InternalMesh != NULL)
      {
      points = this->MoveInternalMesh(this->InternalMesh, pointArray);
      if (points == NULL)
        {
        pointArray->Delete();
        return 0;
        }
      }
    else
      {
      points = vtkPoints::New();
      points->SetData(pointArray);
      }

    // Check if Internal Surface Mesh exists first....
    if (this->SurfaceMesh != NULL)
      {
      if (!this->MoveSurfaceMesh(this->SurfaceMesh, pointArray))
        {
        points->Delete();
        pointArray->Delete();
        return 0;
        }
      }

    if (this->PointZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->PointZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkPolyData::SafeDownCast(this->PointZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    if (this->FaceZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->FaceZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkPolyData::SafeDownCast(this->FaceZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    if (this->CellZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->CellZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkUnstructuredGrid::SafeDownCast(this->CellZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    points->Delete();
    }

  if (createEulerians && moveBoundaryPoints)
    {
    // Check if Boundary Mesh exists first....
    if (this->BoundaryMesh != NULL)
      {
      this->MoveBoundaryMesh(this->BoundaryMesh, pointArray);
      }
    }

  if (pointArray != NULL)
    {
    pointArray->Delete();
    }
  this->Parent->UpdateProgress(0.5);

  if (updateVariables)
    {
    if (createEulerians)
      {
      if (!recreateInternalMesh)
        {
        // Check if Internal Mesh Exists first...
        if (this->InternalMesh != NULL)
          {
          // clean up arrays of the previous timestep
          this->InternalMesh->GetCellData()->Initialize();
          this->InternalMesh->GetPointData()->Initialize();
          }
        if (this->SurfaceMesh != NULL)
          {
          this->SurfaceMesh->GetCellData()->Initialize();
          }
        }
      if (!recreateBoundaryMesh && this->BoundaryMesh != NULL)
        {
        for (unsigned int i = 0; i < this->BoundaryMesh->GetNumberOfBlocks(); i++)
          {
          vtkPolyData *bm =
              vtkPolyData::SafeDownCast(this->BoundaryMesh->GetBlock(i));
          bm->GetCellData()->Initialize();
          bm->GetPointData()->Initialize();
          }
        }
      // read field data variables into Internal/Boundary meshes
      for (int i = 0; i < (int)this->VolFieldFiles->GetNumberOfValues(); i++)
        {
        this->GetVolFieldAtTimeStep(this->InternalMesh, this->BoundaryMesh,
            this->VolFieldFiles->GetValue(i));
        this->Parent->UpdateProgress(0.5 + 0.25 * ((float)(i + 1)
            / ((float)this->VolFieldFiles->GetNumberOfValues() + 0.0001)));
        }
      if (this->SurfaceMesh != NULL || (this->BoundaryMesh != NULL
          && this->BoundaryMesh->GetNumberOfBlocks() > 0))
        {
        for (int i = 0; i < (int)this->SurfaceFieldFiles->GetNumberOfValues();
             i++)
          {
          this->GetSurfaceFieldAtTimeStep(this->SurfaceMesh,
              this->BoundaryMesh, this->SurfaceFieldFiles->GetValue(i));
          this->Parent->UpdateProgress(0.75 + 0.125 * ((float)(i + 1) /
                  ((float)this->SurfaceFieldFiles->GetNumberOfValues()
                  + 0.0001)));
          }
        }
      for (int i = 0; i < (int)this->PointFieldFiles->GetNumberOfValues(); i++)
        {
        this->GetPointFieldAtTimeStep(this->InternalMesh, this->BoundaryMesh,
            this->PointFieldFiles->GetValue(i));
        this->Parent->UpdateProgress(0.75 + 0.125 * ((float)(i + 1)
            / ((float)this->PointFieldFiles->GetNumberOfValues() + 0.0001)));
        }
      }
    }

  if (recreateLagrangianMesh)
    {
    // read lagrangian mesh and fields
    this->LagrangianMesh = this->MakeLagrangianMesh();
    }

  // Add Internal Mesh to final output only if selected for display
  if (this->InternalMesh != NULL)
    {
    output->SetBlock(0, this->InternalMesh);
    this->SetBlockName(output, 0,
        vtkOFFReaderPrivate::InternalMeshIdentifier);
    }

  // Add Internal Surface Mesh to final output only if selected for display
  if (this->SurfaceMesh != NULL)
    {
    const unsigned int groupTypeI = output->GetNumberOfBlocks();
    output->SetBlock(groupTypeI, this->SurfaceMesh);
    this->SetBlockName(output, groupTypeI,
        vtkOFFReaderPrivate::SurfaceMeshIdentifier);
    }

  // set boundary meshes/data as output
  if (this->BoundaryMesh != NULL && this->BoundaryMesh->GetNumberOfBlocks() > 0)
    {
    const unsigned int groupTypeI = output->GetNumberOfBlocks();
    output->SetBlock(groupTypeI, this->BoundaryMesh);
    this->SetBlockName(output, groupTypeI, "[Patches]");
    }

  // set lagrangian mesh as output
  if (this->LagrangianMesh && this->LagrangianMesh->GetNumberOfBlocks() > 0)
    {
    const unsigned int groupTypeI = output->GetNumberOfBlocks();
    output->SetBlock(groupTypeI, this->LagrangianMesh);
    this->SetBlockName(output, groupTypeI, "[Lagrangian Particles]");
    }

  if (this->Parent->GetReadZones())
    {
    vtkMultiBlockDataSet *zones = NULL;
    // set Zone Meshes as output
    if (this->PointZoneMesh != NULL)
      {
      zones = vtkMultiBlockDataSet::New();
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->PointZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "[pointZones]");
      }

    if (this->FaceZoneMesh != NULL)
      {
      if (zones == NULL)
        {
        zones = vtkMultiBlockDataSet::New();
        }
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->FaceZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "[faceZones]");
      }

    if (this->CellZoneMesh != NULL)
      {
      if (zones == NULL)
        {
        zones = vtkMultiBlockDataSet::New();
        }
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->CellZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "[cellZones]");
      }
    if (zones != NULL)
      {
      const unsigned int groupTypeI = output->GetNumberOfBlocks();
      output->SetBlock(groupTypeI, zones);
      this->SetBlockName(output, groupTypeI, "[Zones]");
      }
    }

  if (this->Parent->GetCacheMesh())
    {
    this->TimeStepOld = this->TimeStep;
    }
  else
    {
    this->ClearMeshes();
    this->TimeStepOld = -1;
    }
  this->InternalMeshSelectionStatusOld = this->InternalMeshSelectionStatus;
  this->SurfaceMeshSelectionStatusOld = this->SurfaceMeshSelectionStatus;

  this->Parent->UpdateProgress(1.0);
  return 1;
}

//-----------------------------------------------------------------------------
// static member functions
int vtkOFFReader::GetBoundaryTypeProcessor()
{
  return vtkOFFReaderPrivate::GetBoundaryTypeProcessor();
}

int vtkOFFReader::GetBoundaryTypeInternal()
{
  return vtkOFFReaderPrivate::GetBoundaryTypeInternal();
}

//-----------------------------------------------------------------------------
// constructor
vtkOFFReader::vtkOFFReader()
{
  this->SetNumberOfInputPorts(0);

  this->Parent = this;
  // must be false to avoid reloading by vtkAppendCompositeDataLeaves::Update()
  this->Refresh = false;

  // INTIALIZE FILE NAME
  this->FileName = NULL;
  this->FileNameOld = new vtkStdString;

  // Case path
  this->CasePath = vtkCharArray::New();

  // Child readers
  this->Readers = vtkCollection::New();

  // VTK CLASSES
  this->PatchDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->SurfaceDataArraySelection = vtkDataArraySelection::New();
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->LagrangianDataArraySelection = vtkDataArraySelection::New();

  this->PatchSelectionMTimeOld = 0;
  this->CellSelectionMTimeOld = 0;
  this->SurfaceSelectionMTimeOld = 0;
  this->PointSelectionMTimeOld = 0;
  this->LagrangianSelectionMTimeOld = 0;
  this->LagrangianPathsMTimeOld = 0;

  // for creating cell-to-point translated data
  this->CreateCellToPoint = 1;
  this->CreateCellToPointOld = 1;

  // for caching mesh
  this->CacheMesh = 1;

  // for decomposing polyhedra
  this->DecomposePolyhedra = 1;
  this->DecomposePolyhedraOld = 1;

  // for reading old binary lagrangian/positions format
  this->PositionsIsIn13Format = 0; // turned off by default
  this->PositionsIsIn13FormatOld = 0;

  // for reading single precision binary format
  this->IsSinglePrecisionBinary = 0; // turned off by default
  this->IsSinglePrecisionBinaryOld = 0;

  // for reading zones
  this->ReadZones = 0; // turned off by default
  this->ReadZonesOld = 0;

  // whether to output processor patches or not
  this->OutputProcessorPatches = PROCESSOR_PATCHES_ON;
  this->OutputProcessorPatchesOld = PROCESSOR_PATCHES_ON;

  // determine if time directories are to be listed according to controlDict
  this->ListTimeStepsByControlDict = 0;
  this->ListTimeStepsByControlDictOld = 0;

  // add dimensions to array names
  this->AddDimensionsToArrayNames = 0;
  this->AddDimensionsToArrayNamesOld = 0;

  // has mesh changed at this time step or not
  this->MeshChanged = 0;

  // Lagrangian paths
  this->LagrangianPaths = vtkStringArray::New();

  this->CurrentReaderIndex = 0;
  this->NumberOfReaders = 0;
}

//-----------------------------------------------------------------------------
// destructor
vtkOFFReader::~vtkOFFReader()
{
  this->LagrangianPaths->Delete();

  this->PatchDataArraySelection->Delete();
  this->CellDataArraySelection->Delete();
  this->SurfaceDataArraySelection->Delete();
  this->PointDataArraySelection->Delete();
  this->LagrangianDataArraySelection->Delete();

  this->Readers->Delete();
  this->CasePath->Delete();

  this->SetFileName(0);
  delete this->FileNameOld;
}

//-----------------------------------------------------------------------------
// CanReadFile
int vtkOFFReader::CanReadFile(const char *vtkNotUsed(fileName))
{
  return 1; // so far CanReadFile does nothing.
}

//-----------------------------------------------------------------------------
// PrintSelf
void vtkOFFReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
      << endl;
  os << indent << "Refresh: " << this->Refresh << endl;
  os << indent << "CreateCellToPoint: " << this->CreateCellToPoint << endl;
  os << indent << "CacheMesh: " << this->CacheMesh << endl;
  os << indent << "DecomposePolyhedra: " << this->DecomposePolyhedra << endl;
  os << indent << "PositionsIsIn13Format: " << this->PositionsIsIn13Format
      << endl;
  os << indent << "IsSinglePrecisionBinary: " << this->IsSinglePrecisionBinary
      << endl;
  os << indent << "ReadZones: " << this->ReadZones << endl;
  os << indent << "OutputProcessorPatches: " << this->OutputProcessorPatches << endl;
  os << indent << "ListTimeStepsByControlDict: "
      << this->ListTimeStepsByControlDict << endl;
  os << indent << "AddDimensionsToArrayNames: "
      << this->AddDimensionsToArrayNames << endl;

  os << indent << "Case Path: \n";
  this->CasePath->PrintSelf(os, indent.GetNextIndent());

  this->Readers->InitTraversal();
  vtkObject *reader;
  while ((reader = this->Readers->GetNextItemAsObject()) != NULL)
    {
    os << indent << "Reader instance " << static_cast<void *>(reader) << ": \n";
    reader->PrintSelf(os, indent.GetNextIndent());
    }

  os << indent << "Patch Data Array Selection: \n";
  this->PatchDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Cell Data Array Selection: \n";
  this->CellDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Surface Data Array Selection: \n";
  this->SurfaceDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Point Data Array Selection: \n";
  this->PointDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Lagrangian Data Array Selection: \n";
  this->LagrangianDataArraySelection->PrintSelf(os, indent.GetNextIndent());

  os << indent << "Patch Selection MTime Old: "
    << this->PatchSelectionMTimeOld << endl;
  os << indent << "Cell Selection MTime Old: "
    << this->CellSelectionMTimeOld << endl;
  os << indent << "Surface Selection MTime Old: "
    << this->SurfaceSelectionMTimeOld << endl;
  os << indent << "Point Selection MTime Old: "
    << this->PointSelectionMTimeOld << endl;
  os << indent << "Lagrangian Selection MTime Old: "
    << this->LagrangianSelectionMTimeOld << endl;
  os << indent << "Lagrangian Paths MTime Old: "
    << this->LagrangianPathsMTimeOld << endl;

  os << indent << "FileNameOld: " << (this->FileNameOld ? *this->FileNameOld :
      vtkStdString()) << endl;
  os << indent << "CreateCellToPointOld: " << this->CreateCellToPointOld
      << endl;
  os << indent << "DecomposePolyhedraOld: " << this->DecomposePolyhedraOld
      << endl;
  os << indent << "PositionsIsIn13FormatOld: " << this->PositionsIsIn13FormatOld
      << endl;
  os << indent << "IsSinglePrecisionBinaryOld: "
      << this->IsSinglePrecisionBinaryOld << endl;
  os << indent << "ReadZonesOld: " << this->ReadZonesOld << endl;
  os << indent << "OutputProcessorPatchesOld: " << this->OutputProcessorPatchesOld << endl;
  os << indent << "ListTimeStepsByControlDictOld: "
      << this->ListTimeStepsByControlDictOld << endl;
  os << indent << "AddDimensionsToArrayNamesOld: "
      << this->AddDimensionsToArrayNamesOld << endl;
  os << indent << "MeshChanged: " << this->MeshChanged << endl;

  os << indent << "Lagrangian Paths: \n";
  this->LagrangianPaths->PrintSelf(os, indent.GetNextIndent());

  os << indent << "Number Of Readers: " << this->NumberOfReaders << endl;
  os << indent << "Current Reader Index: " << this->CurrentReaderIndex << endl;

  os << indent << "Parent reader instance: "
      << static_cast<void *>(this->Parent) << endl;
  return;
}

//-----------------------------------------------------------------------------
// selection list handlers

int vtkOFFReader::GetNumberOfSelectionArrays(vtkDataArraySelection *s)
{
  return s->GetNumberOfArrays();
}

int vtkOFFReader::GetSelectionArrayStatus(vtkDataArraySelection *s,
    const char *name)
{
  return s->ArrayIsEnabled(name);
}

void vtkOFFReader::SetSelectionArrayStatus(vtkDataArraySelection *s,
    const char* name, int status)
{
  unsigned long int mTime = s->GetMTime();
  if (status)
    {
    s->EnableArray(name);
    }
  else
    {
    s->DisableArray(name);
    }
  if (mTime != s->GetMTime()) // indicate that the pipeline needs to be updated
    {
    this->Modified();
    }
}

const char *vtkOFFReader::GetSelectionArrayName(vtkDataArraySelection *s,
    int index)
{
  return s->GetArrayName(index);
}

void vtkOFFReader::DisableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->DisableAllArrays();
  if (mTime != s->GetMTime())
    {
    this->Modified();
    }
}

void vtkOFFReader::EnableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->EnableAllArrays();
  if (mTime != s->GetMTime())
    {
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
// RequestInformation
int vtkOFFReader::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
  this->Parent->SetErrorCode(vtkErrorCode::NoError);

  if (!this->FileName || strlen(this->FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    this->Parent->SetErrorCode(vtkErrorCode::NoFileNameError);
    return 0;
    }

  if (this->Parent == this && (*this->FileNameOld != this->FileName
      || this->ListTimeStepsByControlDict
          != this->ListTimeStepsByControlDictOld || this->Refresh))
    {
    // retain selection status when just refreshing a case
    if (*this->FileNameOld != "" && *this->FileNameOld != this->FileName)
      {
      // clear selections
      this->CellDataArraySelection->RemoveAllArrays();
      this->SurfaceDataArraySelection->RemoveAllArrays();
      this->PointDataArraySelection->RemoveAllArrays();
      this->LagrangianDataArraySelection->RemoveAllArrays();
      this->PatchDataArraySelection->RemoveAllArrays();
      }

    // Reset NumberOfReaders here so that the variable will not be
    // reset unwantedly when MakeInformationVector() is called from
    // vtkPOpenFOAMReader
    this->NumberOfReaders = 0;

    if (!this->MakeInformationVector(outputVector, vtkStdString(""))
        || !this->MakeMetaDataAtTimeStep(true))
      {
      return 0;
      }
    this->Refresh = false;
    }
  return 1;
}

//-----------------------------------------------------------------------------
// RequestData
int vtkOFFReader::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
  this->Parent->SetErrorCode(vtkErrorCode::NoError);

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet
      *output =
          vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int nSteps = 0;
  double *requestedTimeValues = NULL;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    requestedTimeValues
        = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }

  if (nSteps > 0)
    {
    outInfo->Set(vtkDataObject::DATA_TIME_STEPS(), requestedTimeValues, 1);
    this->SetTimeValue(requestedTimeValues[0]);
    }

  if (this->Parent == this)
    {
    output->GetFieldData()->AddArray(this->CasePath);
    if (!this->MakeMetaDataAtTimeStep(false))
      {
      return 0;
      }
    this->MeshChanged = 0;
    this->CurrentReaderIndex = 0;
    }

  // compute flags
  // internal mesh selection change is detected within each reader
  const bool recreateInternalMesh = (!this->Parent->CacheMesh)
      || this->Parent->DecomposePolyhedra
          != this->Parent->DecomposePolyhedraOld || this->Parent->ReadZones
      != this->Parent->ReadZonesOld || this->Parent->ListTimeStepsByControlDict
      != this->Parent->ListTimeStepsByControlDictOld
      || this->Parent->IsSinglePrecisionBinary
      != this->Parent->IsSinglePrecisionBinaryOld;
  const bool recreateBoundaryMesh =
      this->Parent->PatchDataArraySelection->GetMTime()
          != this->Parent->PatchSelectionMTimeOld
          || this->Parent->CreateCellToPoint
              != this->Parent->CreateCellToPointOld
          || this->Parent->OutputProcessorPatches
              != this->Parent->OutputProcessorPatchesOld;
  const bool updateVariables = this->Parent->CellDataArraySelection->GetMTime()
      != this->Parent->CellSelectionMTimeOld
      || this->Parent->SurfaceDataArraySelection->GetMTime()
          != this->Parent->SurfaceSelectionMTimeOld
      || this->Parent->PointDataArraySelection->GetMTime()
          != this->Parent->PointSelectionMTimeOld
      || this->Parent->AddDimensionsToArrayNames
          != this->Parent->AddDimensionsToArrayNamesOld;
  const bool recreateLagrangianMesh = (!this->Parent->CacheMesh)
      || this->Parent->PatchDataArraySelection->GetMTime()
          != this->Parent->PatchSelectionMTimeOld
      || this->Parent->LagrangianDataArraySelection->GetMTime()
          != this->Parent->LagrangianSelectionMTimeOld
      || this->Parent->LagrangianPaths->GetMTime()
          != this->Parent->LagrangianPathsMTimeOld
      || this->Parent->ListTimeStepsByControlDict
          != this->Parent->ListTimeStepsByControlDictOld
      || this->Parent->CreateCellToPoint
          != this->Parent->CreateCellToPointOld
      || this->Parent->PositionsIsIn13Format
          != this->Parent->PositionsIsIn13FormatOld
      || this->Parent->IsSinglePrecisionBinary
          != this->Parent->IsSinglePrecisionBinaryOld;

  // create dataset
  int ret = 1;
  vtkOFFReaderPrivate *reader;
  // if the only region is not a subregion, omit being wrapped by a
  // multiblock dataset
  if (this->Readers->GetNumberOfItems() == 1 && (reader = vtkOFFReaderPrivate::SafeDownCast(
          this->Readers->GetItemAsObject(0)))->GetRegionName() == "")
    {
    ret = reader->RequestData(output, recreateInternalMesh,
        recreateBoundaryMesh, updateVariables, recreateLagrangianMesh);
    this->Parent->CurrentReaderIndex++;
    }
  else
    {
    this->Readers->InitTraversal();
    while ((reader
        = vtkOFFReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
        != NULL)
      {
      vtkMultiBlockDataSet *subOutput = vtkMultiBlockDataSet::New();
      if (reader->RequestData(subOutput, recreateInternalMesh,
          recreateBoundaryMesh, updateVariables, recreateLagrangianMesh))
        {
        vtkStdString regionName(reader->GetRegionName());
        if (regionName == "")
          {
          regionName = "defaultRegion";
          }
        const int blockI = output->GetNumberOfBlocks();
        output->SetBlock(blockI, subOutput);
        output->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), regionName.c_str());
        }
      else
        {
        ret = 0;
        }
      subOutput->Delete();
      this->Parent->CurrentReaderIndex++;
      }
    }

  if (this->Parent == this) // update only if this is the top-level reader
    {
    this->UpdateStatus();
    }

  return ret;
}

//-----------------------------------------------------------------------------
void vtkOFFReader::SetTimeInformation(vtkInformationVector *outputVector,
    vtkDoubleArray *timeValues)
{
  // the length of the time value information array should not be zero
  // -- otherwise ParaView will crash.
  if (timeValues->GetNumberOfTuples() > 0)
    {
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        timeValues->GetPointer(0), timeValues->GetNumberOfTuples());

    double timeRange[2];
    timeRange[0] = timeValues->GetValue(0);
    timeRange[1] = timeValues->GetValue(timeValues->GetNumberOfTuples() - 1);
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReader::GetRegions(vtkStringArray *regionNames,
    const vtkStdString &timeDir)
{
  vtkDirectory *dir = vtkDirectory::New();
  if (!dir->Open(timeDir.c_str()))
    {
    dir->Delete();
    return;
    }
  for (int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
    {
    vtkStdString subDir(dir->GetFile(fileI));
    // "uniform" subdir is used for recording time by OF <= 1.4.1
    if (dir->FileIsDirectory(subDir.c_str()) && subDir != "." && subDir != ".."
        && subDir != "polyMesh" && subDir != "lagrangian" && subDir != "uniform")
      {
      if(regionNames->LookupValue(subDir) >= 0)
        {
        continue;
        }
      vtkStdString boundaryPath(timeDir + "/" + subDir + "/polyMesh/boundary");
      if (vtksys::SystemTools::FileExists(boundaryPath.c_str(), true)
          || vtksys::SystemTools::FileExists((boundaryPath + ".gz").c_str(), true))
        {
        regionNames->InsertNextValue(subDir);
        }
      }
    }
  dir->Delete();
}

//-----------------------------------------------------------------------------
int vtkOFFReader::MakeInformationVector(
    vtkInformationVector *outputVector, const vtkStdString& procName)
{
  *this->FileNameOld = vtkStdString(this->FileName);

  // clear prior case information
  this->Readers->RemoveAllItems();

  // recreate case information
  vtkStdString casePath, controlDictPath;
  this->CreateCasePath(casePath, controlDictPath);
  casePath += procName + (procName == "" ? "" : "/");
  vtkOFFReaderPrivate *masterReader = vtkOFFReaderPrivate::New();
  if (!masterReader->MakeInformationVector(casePath, controlDictPath, procName,
      this->Parent))
    {
    masterReader->Delete();
    return 0;
    }

  if (outputVector != NULL)
    {
    this->SetTimeInformation(outputVector, masterReader->GetTimeValues());
    }

  if (masterReader->GetTimeValues()->GetNumberOfTuples() == 0)
    {
    vtkWarningMacro(<< "Case " << casePath.c_str()
        << " contains no timestep data.");
    masterReader->Delete();
    this->Parent->SetErrorCode(vtkErrorCode::NoError);
    return 1;
    }

  this->Readers->AddItem(masterReader);

  // search subregions under constant subdirectory
  vtkStringArray *regionNames = vtkStringArray::New();
  this->GetRegions(regionNames, casePath + "constant");
  vtkStringArray *timeNames = masterReader->GetTimeNames();
  int nTimes = timeNames->GetNumberOfValues() > 2
      ? 2 : timeNames->GetNumberOfValues();
  for(int timeI = 0; timeI < nTimes; timeI++)
    {
    this->GetRegions(regionNames, casePath + timeNames->GetValue(timeI));
    }
  regionNames->Squeeze();
  for(int regionI = 0; regionI < regionNames->GetNumberOfValues(); regionI++)
    {
    vtkOFFReaderPrivate *subReader = vtkOFFReaderPrivate::New();
    subReader->SetupInformation(casePath, regionNames->GetValue(regionI),
        procName, masterReader);
    this->Readers->AddItem(subReader);
    subReader->Delete();
    }
  regionNames->Delete();
  masterReader->Delete();
  this->Parent->NumberOfReaders += this->Readers->GetNumberOfItems();

  if (this->Parent == this)
    {
    this->CreateCharArrayFromString(this->CasePath, "CasePath", casePath);
    }

  return 1;
}

//-----------------------------------------------------------------------------
void vtkOFFReader::CreateCasePath(vtkStdString &casePath,
    vtkStdString &controlDictPath)
{
#if defined(_WIN32)
  const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
  const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
  controlDictPath = this->FileName;

  // determine the case directory and path to controlDict
  vtkStdString::size_type pos = controlDictPath.find_last_of(pathFindSeparator);
  if (pos == vtkStdString::npos)
    {
    // if there's no prepending path, prefix with the current directory
    controlDictPath = "." + pathSeparator + controlDictPath;
    pos = 1;
    }
  // filenames can be specified by wildcards as well as the classic
  // ".foam" extensions.
  // cf. vtkSMReaderFactory::vtkInternals::vtkValue::FillInformation()
  vtkStdString lastCmpt(controlDictPath.substr(pos + 1));
  if (lastCmpt == "controlDict" || lastCmpt == "fvSchemes"
      || lastCmpt == "fvSolution")
    {
    // remove trailing filename
    casePath = controlDictPath.substr(0, pos);
    if (casePath == ".")
      {
      casePath = ".." + pathSeparator;
      }
    else
      {
      vtkStdString::size_type pos2 = casePath.find_last_of(pathFindSeparator);
      if (pos2 == vtkStdString::npos)
        {
        casePath = "." + pathSeparator;
        }
      else
        {
        if (casePath.substr(pos2 + 1) == "system")
          {
          // remove trailing "system"
          casePath.erase(pos2 + 1); // preserve the last "/"
          }
        else
          {
          // maybe a multiregion case; go up one more level blindly
          // assuming the directory name is "system"
          vtkStdString multiRegionCasePath(casePath.substr(0, pos2));
          vtkStdString::size_type pos3
              = multiRegionCasePath.find_last_of(pathFindSeparator);
          if (pos3 != vtkStdString::npos)
            {
            casePath.erase(pos3 + 1);
            }
          }
        }
      }
    }
  else
    {
    // if the file is named other than controlDict / fvSchemes /
    // fvSolution, use the directory containing the file as case
    // directory
    casePath = controlDictPath.substr(0, pos + 1);
    controlDictPath = casePath + "system" + pathSeparator + "controlDict";
    }
}

//-----------------------------------------------------------------------------
void vtkOFFReader::AddSelectionNames(vtkDataArraySelection *selections,
    vtkStringArray *objects)
{
  objects->Squeeze();
  vtkSortDataArray::Sort(objects);
  for (int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->AddArray(objects->GetValue(nameI).c_str());
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
bool vtkOFFReader::SetTimeValue(const double timeValue)
{
  bool modified = false;
  vtkOFFReaderPrivate *reader;
  this->Readers->InitTraversal();
  while ((reader
      = vtkOFFReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
      != NULL)
    {
    const unsigned long mTime = reader->GetMTime();
    reader->SetTimeValue(timeValue);
    if (reader->GetMTime() != mTime)
      {
      modified = true;
      }
    }
  return modified;
}

//-----------------------------------------------------------------------------
vtkDoubleArray *vtkOFFReader::GetTimeValues()
{
  if (this->Readers->GetNumberOfItems() <= 0)
    {
    return NULL;
    }
  vtkOFFReaderPrivate *reader =
      vtkOFFReaderPrivate::SafeDownCast(this->Readers->GetItemAsObject(0));
  return reader != NULL ? reader->GetTimeValues() : NULL;
}

//-----------------------------------------------------------------------------
int vtkOFFReader::MakeMetaDataAtTimeStep(const bool listNextTimeStep)
{
  vtkStringArray *cellSelectionNames = vtkStringArray::New();
  vtkStringArray *surfaceSelectionNames = vtkStringArray::New();
  vtkStringArray *pointSelectionNames = vtkStringArray::New();
  vtkStringArray *lagrangianSelectionNames = vtkStringArray::New();
  int ret = 1;
  vtkOFFReaderPrivate *reader;
  this->Readers->InitTraversal();
  while ((reader
      = vtkOFFReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
      != NULL)
    {
    ret *= reader->MakeMetaDataAtTimeStep(cellSelectionNames,
        surfaceSelectionNames, pointSelectionNames, lagrangianSelectionNames,
        listNextTimeStep);
    }
  this->AddSelectionNames(this->Parent->CellDataArraySelection,
      cellSelectionNames);
  this->AddSelectionNames(this->Parent->SurfaceDataArraySelection,
      surfaceSelectionNames);
  this->AddSelectionNames(this->Parent->PointDataArraySelection,
      pointSelectionNames);
  this->AddSelectionNames(this->Parent->LagrangianDataArraySelection,
      lagrangianSelectionNames);

  return ret;
}

//-----------------------------------------------------------------------------
void vtkOFFReader::CreateCharArrayFromString(vtkCharArray *array,
    const char *name, vtkStdString &string)
{
  array->Initialize();
  array->SetName(name);
  const size_t len = string.length();
  char *ptr = array->WritePointer(0, static_cast<vtkIdType>(len + 1));
  memcpy(ptr, string.c_str(), len);
  ptr[len] = '\0';
}

//-----------------------------------------------------------------------------
void vtkOFFReader::UpdateStatus()
{
  // update selection MTimes
  this->PatchSelectionMTimeOld = this->PatchDataArraySelection->GetMTime();
  this->CellSelectionMTimeOld = this->CellDataArraySelection->GetMTime();
  this->SurfaceSelectionMTimeOld = this->SurfaceDataArraySelection->GetMTime();
  this->PointSelectionMTimeOld = this->PointDataArraySelection->GetMTime();
  this->LagrangianSelectionMTimeOld
      = this->LagrangianDataArraySelection->GetMTime();
  this->LagrangianPathsMTimeOld = this->LagrangianPaths->GetMTime();
  this->CreateCellToPointOld = this->CreateCellToPoint;
  this->DecomposePolyhedraOld = this->DecomposePolyhedra;
  this->PositionsIsIn13FormatOld = this->PositionsIsIn13Format;
  this->IsSinglePrecisionBinaryOld = this->IsSinglePrecisionBinary;
  this->ReadZonesOld = this->ReadZones;
  this->OutputProcessorPatchesOld = this->OutputProcessorPatches;
  this->ListTimeStepsByControlDictOld = this->ListTimeStepsByControlDict;
  this->AddDimensionsToArrayNamesOld = this->AddDimensionsToArrayNames;
}

//-----------------------------------------------------------------------------
void vtkOFFReader::UpdateProgress(double amount)
{
  this->vtkAlgorithm::UpdateProgress((static_cast<double>(this->Parent->CurrentReaderIndex)
      + amount) / static_cast<double>(this->Parent->NumberOfReaders));
}
