pipe_ra = '${units 25. cm -> m}'
pipe_dia = '${fparse pipe_ra * 2}'
A_pipe = '${fparse 0.25 * pi * pipe_dia^2}'

T_amb = 293.15 # K
T_in = 570.15 # K
m_dot_in = 1.33 # kg/s
pressure = 15.5e6 # Pa


p_tank_initial = 1e5


[GlobalParams]
  initial_p = ${pressure}
  initial_vel = 0.0001
  initial_T = ${T_in}

  gravity_vector = '0 0 0'

  rdg_slope_reconstruction = minmod
  scaling_factor_1phase = '1 1e-2 1e-4'

  closures = simple_closures
  fp = water
[]


[FluidProperties]
  [water]
    type = SimpleFluidProperties
    molar_mass = 1.8e-2                          # Molar mass of water (kg/mol)
    thermal_expansion = 2.14e-4                  # Coefficient of thermal expansion (1/K)
    cp = 4194.0                                  # Specific heat capacity at constant pressure (J/kg/K)
    cv = 4186.0                                  # Specific heat capacity at constant volume (J/kg/K)
    bulk_modulus = 2.2e9                         # Bulk modulus of water (Pa)
    thermal_conductivity = 0.6                   # Thermal conductivity (W/m/K)
    viscosity = 1.0e-3                           # Dynamic viscosity of water (Pa.s)
    density0 = 1000.0                            # Reference density of water (kg/m^3)
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseTHM
  []
[]

[FunctorMaterials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'density_functor mu'
    prop_values = '1  1'
  []
[]

[AuxKernels]
  [speed]
    type = ParsedAux
    variable = 'velocity_norm'
    coupled_variables = 'superficial_vel_x superficial_vel_y porosity'
    expression = 'sqrt(superficial_vel_x*superficial_vel_x + superficial_vel_y*superficial_vel_y) / porosity'
  []
[]

[AuxVariables]
  [porosity]
    type = MooseVariableFVReal
    initial_condition = 0.5
  []
  [velocity_norm]
    type = MooseVariableFVReal
    block = 'iwst:iwst'
  []
[]

[Physics]
  [NavierStokes]
    [Flow]
      [flow]
        compressibility = 'incompressible'
        porous_medium_treatment = true

        density = 'density_functor'
        dynamic_viscosity = 'mu'

        initial_velocity = '0 0 0'
        initial_pressure = '${p_tank_initial}'

        mass_advection_interpolation = 'upwind'
        momentum_advection_interpolation = 'upwind'

        inlet_boundaries = 'iwst:left'
        momentum_inlet_types = 'fixed-velocity'
        momentum_inlet_function = '0 0 0'

        wall_boundaries = 'iwst:top iwst:bottom iwst:back iwst:front'
        momentum_wall_types = 'noslip symmetry noslip symmetry'

        outlet_boundaries = 'iwst:right'
        momentum_outlet_types = 'fixed-pressure'
        pressure_function = '${p_tank_initial}'
      []
    []
  []
[]



[Components]

  #IWST tank which has been created in a separated directory, meshed directly using the command: orca IWST.i --mesh-only
  [iwst]
    type = FileMeshPhysicsComponent
    file = /root/projects/ORCA/models/otc_simulation/otc_geometries/IWST_in.e
    position = '0 0 0'
    physics = 'flow'
  []

  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'iwst_pipe_1:in'
    m_dot = ${m_dot_in}
    T = ${T_in}
  []
  [iwst_pipe_1]
    type = FlowChannel1Phase
    position = '2.5 2.5 -5'
    orientation = '0 0 1'
    length = 20
    n_elems = 20
    A = ${A_pipe}
    D_h = ${pipe_dia}
  []
  [to_hs]
    type = HeatTransferFromHeatStructure1Phase
    flow_channel = iwst_pipe_1
    hs = hs
    hs_side = inner
    P_hf = '${fparse pi * 2 * pipe_ra * 19}'
  []
  [hs]
    type = HeatStructureCylindrical
    position = '2.5 2.5 -5'
    orientation = '0 0 1'
    length = 20
    n_elems = 20

    inner_radius = ${pipe_ra}
    num_rods = 19
    initial_T = ${T_amb}
    names = 'hs_steel'
    widths = 0.5
    materials = 'steel'
    n_part_elems = '20'
    offset_mesh_by_inner_radius = true
  []
  [from_hs]
    type = HeatTransferFromHeatStructure1Phase
    flow_channel = iwst
    hs = hs
    hs_side = outer
    P_hf = '${fparse pi * 2 * pipe_ra * 19}'
  []

  [outlet]
    type = Outlet1Phase
    input = 'iwst_pipe_1:out'
    p = ${pressure}
  []
[]

[HeatStructureMaterials]
  [steel]
    type = SolidMaterialProperties
    rho = 8050 # kg/m3
    k = 45 # W/(m.K)
    cp = 466 # J/(kg.K)
  []
[]


[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  end_time = 160000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1000.
    growth_factor = 1.2
    cutback_factor = 0.8
  []
  
  dtmin = 1e-5

  solve_type = NEWTON
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7
  nl_max_its = 25

  l_tol = 1e-6
  l_max_its = 10

  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu     '
  line_search = basic
[]


[Outputs]
  exodus = true
  file_base = /root/projects/ORCA/models/otc_simulation/otc_output_files/otc_output
  [console]
    type = Console
    max_rows = 1
    outlier_variable_norms = false
  []
  print_linear_residuals = false
[]