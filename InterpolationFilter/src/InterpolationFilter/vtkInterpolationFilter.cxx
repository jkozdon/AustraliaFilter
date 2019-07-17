#include "vtkInterpolationFilter.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataObject.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLagrangeHexahedron.h"
#include "vtkLagrangeQuadrilateral.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

vtkStandardNewMacro(vtkInterpolationFilter);

int vtkInterpolationFilter::number_of_row_points(vtkCell *cell)
{
  int Np = cell->GetNumberOfPoints();
  int Nq = 0;
  switch (cell->GetCellType())
  {
  case VTK_LAGRANGE_HEXAHEDRON:
  {
    switch (Np)
    {
    case 2 * 2 * 2:
    {
      Nq = 2;
      break;
    }
    case 3 * 3 * 3:
    {
      Nq = 3;
      break;
    }
    case 4 * 4 * 4:
    {
      Nq = 4;
      break;
    }
    case 5 * 5 * 5:
    {
      Nq = 5;
      break;
    }
    case 6 * 6 * 6:
    {
      Nq = 6;
      break;
    }
    case 7 * 7 * 7:
    {
      Nq = 7;
      break;
    }
    case 8 * 8 * 8:
    {
      Nq = 8;
      break;
    }
    case 9 * 9 * 9:
    {
      Nq = 9;
      break;
    }
    case 10 * 10 * 10:
    {
      Nq = 10;
      break;
    }
    }
    if (Nq * Nq * Nq != Np)
    {
      vtkErrorMacro("problem determine order for VTK_LAGRANGE_HEXAHEDRON");
    }
    break;
  }

  case VTK_LAGRANGE_QUADRILATERAL:
  {
    {
      switch (Np)
      {
      case 2 * 2:
      {
        Nq = 2;
        break;
      }
      case 3 * 3:
      {
        Nq = 3;
        break;
      }
      case 4 * 4:
      {
        Nq = 4;
        break;
      }
      case 5 * 5:
      {
        Nq = 5;
        break;
      }
      case 6 * 6:
      {
        Nq = 6;
        break;
      }
      case 7 * 7:
      {
        Nq = 7;
        break;
      }
      case 8 * 8:
      {
        Nq = 8;
        break;
      }
      case 9 * 9:
      {
        Nq = 9;
        break;
      }
      case 10 * 10:
      {
        Nq = 10;
        break;
      }
      }
      if (Nq * Nq != Np)
      {
        vtkErrorMacro("problem determine order for VTK_LAGRANGE_QUADRILATERAL");
      }
      break;
    }
    break;
  }
  default:
  {
    vtkErrorMacro("only VTK_LAGRANGE_QUADRILATERAL and "
                  "VTK_LAGRANGE_HEXAHEDRON are supported");
  }
  }
  return Nq;
}

vtkInterpolationFilter::vtkInterpolationFilter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
};
vtkInterpolationFilter::~vtkInterpolationFilter(){};

int vtkInterpolationFilter::RequestData(vtkInformation *vtkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector)
{
  vtkUnstructuredGrid *input = this->GetUnstructuredGridInput(0);
  vtkUnstructuredGrid *output = this->GetOutput();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

#ifndef NDEBUG
  this->DebugOn();
#endif

#define MAX_NUM_COMP 10
  double fld[MAX_NUM_COMP * 10 * 10 * 10];

  output->SetCells(input->GetCellTypesArray(), input->GetCellLocationsArray(),
                   input->GetCells());
  vtkPoints *in_points = input->GetPoints();
  points->DeepCopy(in_points);
  output->SetPoints(points);
  vtkPoints *out_points = output->GetPoints();

#define NUM_DIMENSIONS 3
  for (int e = 0; e < input->GetNumberOfCells(); ++e)
  {
    // Get the cell we are working on
    vtkCell *cell = input->GetCell(e);

    // Determine how many points are in the cell
    int Np = cell->GetNumberOfPoints();

    // Get the list of IDs for this cell
    vtkIdList *ids = cell->GetPointIds();

    // Check that the number of points and IDs are the same
    assert(Np == ids->GetNumberOfIds());

    // Number of row points
    int Nq = number_of_row_points(cell);
    int order[] = {Nq - 1, Nq - 1, Nq - 1};

    if (Nq == 0)
    {
      vtkErrorMacro("only something went wrong in number_of_row_points");
      return 0;
    }

    // Load the data
    switch (cell->GetCellType())
    {
    case VTK_LAGRANGE_HEXAHEDRON:
    {
      for (int k = 0; k < Nq; ++k)
        for (int j = 0; j < Nq; ++j)
          for (int i = 0; i < Nq; ++i)
          {
            int m = i + Nq * (j + k * Nq);
            int n = vtkLagrangeHexahedron::PointIndexFromIJK(i, j, k, order);
            in_points->GetPoint(ids->GetId(n), fld + NUM_DIMENSIONS * m);
          }
      break;
    }
    case VTK_LAGRANGE_QUADRILATERAL:
    {
      for (int j = 0; j < Nq; ++j)
        for (int i = 0; i < Nq; ++i)
        {
          int m = i + Nq * j;
          int n = vtkLagrangeQuadrilateral::PointIndexFromIJK(i, j, order);
          in_points->GetPoint(ids->GetId(n), fld + NUM_DIMENSIONS * m);
        }
      break;
    }
    default:
    {
      vtkErrorMacro("only VTK_LAGRANGE_QUADRILATERAL and "
                    "VTK_LAGRANGE_HEXAHEDRON are supported");
      return 0;
    }
    }

    // Do Interpolation!
    // TODO: Right now just flipping y
    for (int n = 0; n < Np; ++n)
      fld[1 + NUM_DIMENSIONS * n] = -fld[1 + NUM_DIMENSIONS * n];

    // Write out the data
    switch (cell->GetCellType())
    {
    case VTK_LAGRANGE_HEXAHEDRON:
    {
      for (int k = 0; k < Nq; ++k)
        for (int j = 0; j < Nq; ++j)
          for (int i = 0; i < Nq; ++i)
          {
            int m = i + Nq * (j + k * Nq);
            int n = vtkLagrangeHexahedron::PointIndexFromIJK(i, j, k, order);
            out_points->SetPoint(ids->GetId(n), fld + NUM_DIMENSIONS * m);
          }
      break;
    }
    case VTK_LAGRANGE_QUADRILATERAL:
    {
      for (int j = 0; j < Nq; ++j)
        for (int i = 0; i < Nq; ++i)
        {
          int m = i + Nq * j;
          int n = vtkLagrangeQuadrilateral::PointIndexFromIJK(i, j, order);
          out_points->SetPoint(ids->GetId(n), fld + NUM_DIMENSIONS * m);
        }
      break;
    }
    default:
    {
      vtkErrorMacro("only VTK_LAGRANGE_QUADRILATERAL and "
                    "VTK_LAGRANGE_HEXAHEDRON are supported");
      return 0;
    }
    }
  }

  output->GetCellData()->DeepCopy(input->GetCellData());
  output->GetPointData()->DeepCopy(input->GetPointData());

  // For each field we loop through the elements and interpolate
  for (vtkIdType arr = 0; arr < input->GetPointData()->GetNumberOfArrays();
       ++arr)
  {
    // Get the in data array
    vtkDataArray *in_data = input->GetPointData()->GetArray(arr);
    // Get the out data array
    vtkDataArray *out_data = output->GetPointData()->GetArray(arr);
    // Get the number of components in the array data array
    int num_comp = in_data->GetNumberOfComponents();
    assert(num_comp <= MAX_NUM_COMP);
    // Loop through all the elements
    for (int e = 0; e < input->GetNumberOfCells(); ++e)
    {
      // Get the cell we are working on
      vtkCell *cell = input->GetCell(e);

      // Determine how many points are in the cell
      int Np = cell->GetNumberOfPoints();

      // Get the list of IDs for this cell
      vtkIdList *ids = cell->GetPointIds();

      // Check that the number of points and IDs are the same
      assert(Np == ids->GetNumberOfIds());

      // Number of row points
      int Nq = number_of_row_points(cell);
      int order[] = {Nq - 1, Nq - 1, Nq - 1};

      if (Nq == 0)
      {
        vtkErrorMacro("only something went wrong in number_of_row_points");
        return 0;
      }

      // Load the data
      switch (cell->GetCellType())
      {
      case VTK_LAGRANGE_HEXAHEDRON:
      {
        for (int k = 0; k < Nq; ++k)
          for (int j = 0; j < Nq; ++j)
            for (int i = 0; i < Nq; ++i)
            {
              int m = i + Nq * (j + k * Nq);
              int n = vtkLagrangeHexahedron::PointIndexFromIJK(i, j, k, order);
              in_data->GetTuple(ids->GetId(n), fld + num_comp * m);
            }
        break;
      }
      case VTK_LAGRANGE_QUADRILATERAL:
      {
        for (int j = 0; j < Nq; ++j)
          for (int i = 0; i < Nq; ++i)
          {
            int m = i + Nq * j;
            int n = vtkLagrangeQuadrilateral::PointIndexFromIJK(i, j, order);
            in_data->GetTuple(ids->GetId(n), fld + num_comp * m);
          }
        break;
      }
      default:
      {
        vtkErrorMacro("only VTK_LAGRANGE_QUADRILATERAL and "
                      "VTK_LAGRANGE_HEXAHEDRON are supported");
        return 0;
      }
      }

      // Do Interpolation!
      // TODO: Right now just averaging
      for (int c = 0; c < num_comp; ++c)
      {
        double avg = 0;
        for (int n = 0; n < Np; ++n)
          avg += fld[c + n * num_comp];
        for (int n = 0; n < Np; ++n)
          fld[c + n * num_comp] = avg / Np;
      }

      // Write out the data
      switch (cell->GetCellType())
      {
      case VTK_LAGRANGE_HEXAHEDRON:
      {
        for (int k = 0; k < Nq; ++k)
          for (int j = 0; j < Nq; ++j)
            for (int i = 0; i < Nq; ++i)
            {
              int m = i + Nq * (j + k * Nq);
              int n = vtkLagrangeHexahedron::PointIndexFromIJK(i, j, k, order);
              out_data->SetTuple(ids->GetId(n), fld + num_comp * m);
            }
        break;
      }
      case VTK_LAGRANGE_QUADRILATERAL:
      {
        for (int j = 0; j < Nq; ++j)
          for (int i = 0; i < Nq; ++i)
          {
            int m = i + Nq * j;
            int n = vtkLagrangeQuadrilateral::PointIndexFromIJK(i, j, order);
            out_data->SetTuple(ids->GetId(n), fld + num_comp * m);
          }
        break;
      }
      default:
      {
        vtkErrorMacro("only VTK_LAGRANGE_QUADRILATERAL and "
                      "VTK_LAGRANGE_HEXAHEDRON are supported");
        return 0;
      }
      }
    }
  }

  return 1;
}

int vtkInterpolationFilter::RequestUpdateExtent(
    vtkInformation *request, vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  //  this->DebugOn();

  int piece, numPieces, ghostLevels;

  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces =
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevels = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  if (numPieces > 1)
  {
    ++ghostLevels;
  }

  vtkDebugMacro(<< "Running Update Extent" << piece << numPieces
                << ghostLevels);

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
              numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
              ghostLevels);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
};

int vtkInterpolationFilter::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}
