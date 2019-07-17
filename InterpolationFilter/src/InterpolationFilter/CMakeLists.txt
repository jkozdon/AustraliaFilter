set(classes
  vtkInterpolationFilter)

set(private_headers
  vtkInterpolationFilter.h)

vtk_module_add_module(INTERPOLATIONFILTER
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  InterpolationFilter.xml)
