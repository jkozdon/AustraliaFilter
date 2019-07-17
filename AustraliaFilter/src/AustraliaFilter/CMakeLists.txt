set(classes
  vtkAustraliaFilter)

set(private_headers
  vtkAustraliaFilter.h)

vtk_module_add_module(AUSTRALIAFILTER
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  AustraliaFilter.xml)
