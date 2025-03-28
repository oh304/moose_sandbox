# This input file is used to generate a rectangle mesh for other tests. It
# should be run with "--mesh-only rectangle.e".

[Mesh]
  # [gen]
  #   type = GeneratedMeshGenerator
  #   dim = 2
  #   xmin = 0
  #   xmax = 10
  #   ymin = 0
  #   ymax = 1
  #   nx = 20
  #   ny = 5
  # []
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 5
    zmin = 0
    zmax = 11.2
    nx = 5
    ny = 5
    nz = 15
  []
  [rename_block]
    type = RenameBlockGenerator
    input = gen
    old_block = 0
    new_block = 'iwst'
  []
[]
