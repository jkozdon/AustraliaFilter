#include "vtkCell.h"
#include "vtkUnstructuredGridAlgorithm.h"

class vtkInterpolationFilter : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkInterpolationFilter *New();
  vtkTypeMacro(vtkInterpolationFilter, vtkUnstructuredGridAlgorithm);

protected:
  vtkInterpolationFilter();
  ~vtkInterpolationFilter();

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **InputVector,
                          vtkInformationVector *outputVector) override;
  virtual int RequestUpdateExtent(vtkInformation *request,
                                  vtkInformationVector **InputVector,
                                  vtkInformationVector *outputVector) override;
  virtual int FillInputPortInformation(int port, vtkInformation *info) override;

private:
  int number_of_row_points(vtkCell *cell);
};
