# --- Temperature ---
coolant_initial_temperature = ${units 20 degC -> K}
#inconel_pipe_initial_temperature = ${units 20 degC -> K}
inconel_pipe_initial_temperature = 700
reference_temperature = 1000

# --- Geometry ---
#height_of_tank = 0.1                        # [m] 

# --- Salt Properties (T=1000K) ---

#properties associated with water at 293K
rho = 998.2 # Density [kg/m^3]
drho_dT = -0.07  # Temp. Derivative of Density Eq [kg/m^3/K]
alpha = ${fparse -1/rho * drho_dT}  # Thermal expansion coefficient

cp = 4186     # Specific heat capacity [J/kg/K]
mu = 0.001002     # Dynamic viscosity [Pa.s]
k1 = 0.598        # Thermal conductivity [W/m/K]

# --- Inconel Properties (T=373.15K) ---
rho2 = 8440  # density [kg/m^3]
cp2 = 410.0  # specific heat capacity [J/kg/K] -> RTP: cp2 = 410.0   
k2 = 9.8000  # Thermal conductivity [W/m/K]    -> RTP: k2 = 9.8000

# --- Other Constants ---
g = -9.81 # Accel. due to Gravity [kg/m/s^2]

inlet_temp = ${units 20 degC -> K}
inlet_velocity = 0.2
initial_velocity = 0.2
outlet_pressure = 1e5

# --- Output Configuration ---
output_csv_name = 'coolant_channel_GG'

[Mesh]
    [mesh1]
        type = FileMeshGenerator
        file = ./../../geometry/coolant_channel_GG.msh
    []
[]

[UserObjects]
    [rc]
        type = INSFVRhieChowInterpolator
        block = 'coolant'
        u = u
        v = v
        pressure = pressure
    []
[]

[GlobalParams]
    fv = true
    family = MONOMIAL
    order = CONSTANT
    rhie_chow_user_object = 'rc'
    advected_interp_method = 'upwind'
    velocity_interp_method = 'rc'
[]

[Variables]
    [T]                                     
        type = INSFVEnergyVariable
        block = 'coolant'
        initial_condition = ${coolant_initial_temperature}
    []
    [Ts]                                    
        type = INSFVEnergyVariable
        block = 'inconel_top inconel_bottom'
        initial_condition = ${inconel_pipe_initial_temperature}
    []
    [u]                                    
        type = INSFVVelocityVariable
        cache_cell_gradients = false
        block = 'coolant'
        initial_condition = ${initial_velocity}
    []
    [v]                                 
        type = INSFVVelocityVariable
        cache_cell_gradients = false
        block = 'coolant'
        initial_condition = 1e-9
    []
    [pressure]
        type = INSFVPressureVariable
        block = 'coolant'
    []
    
[]

[AuxVariables]
    [mixing_len]
    []
    [velocity_xy]
        type = MooseVariableFVReal
    []
[]

[FVKernels]
    # --- CONSERVATION OF MASS TERM ---
    [mass]
        type = INSFVMassAdvection
        rho = ${rho}
        variable = pressure
    [] 
    [u_time]
        type = INSFVMomentumTimeDerivative
        variable = u
        rho = ${rho}
        momentum_component = 'x'
    []
    [u_advection]
        type = INSFVMomentumAdvection
        variable = u
        rho = ${rho}
        momentum_component = 'x'
    []
    [u_viscosity]
        type = INSFVMomentumDiffusion
        variable = u
        momentum_component = 'x'
        mu = 'mu'
    []
    [u_pressure]
        type = INSFVMomentumPressure
        variable = u
        momentum_component = 'x'
        pressure = pressure
    []
    [u_buoyancy]
        type = INSFVMomentumBoussinesq
        variable = u
        T_fluid = T
        gravity = '0 ${g} 0'
        ref_temperature = ${reference_temperature}
        rho = ${rho}
        momentum_component = 'x'
        alpha_name = ${alpha}
    []
    [v_buoyancy]
        type = INSFVMomentumBoussinesq
        variable = v
        T_fluid = T
        gravity = '0 ${g} 0'
        ref_temperature = ${reference_temperature}
        rho = ${rho}
        momentum_component = 'y'
        alpha_name = ${alpha}
    []

    [v_time]
        type = INSFVMomentumTimeDerivative
        variable = v
        rho = ${rho}
        momentum_component = 'y'
    []
    [v_advection]
        type = INSFVMomentumAdvection
        variable = v
        rho = ${rho}
        momentum_component = 'y'
    []
    [v_viscosity]
        type = INSFVMomentumDiffusion
        variable = v
        momentum_component = 'y'
        mu = 'mu'
    []
    [v_pressure]
        type = INSFVMomentumPressure
        variable = v
        momentum_component = 'y'
        pressure = pressure
    []

    [temp_time]
        type = INSFVEnergyTimeDerivative
        variable = T
        rho = ${rho}
        dh_dt = 'dh_dt'
    []

    [heat_conduction_in_coolant]
        type = FVDiffusion
        variable = T
        #coeff = 'cp'
        coeff = ${k1}
    []
    [heat_conduction_in_inconel]
        type = FVDiffusion
        variable = Ts
        coeff = ${k2}
    []
    [temp_advection_in_coolant]
        type = INSFVEnergyAdvection
        variable = T
    []
    [time_derivative_in_fs]
        type = FVHeatConductionTimeDerivative
        specific_heat = ${cp}
        density_name = ${rho}
        variable = T
    []
    [time_derivative_in_inconel]
        type = FVHeatConductionTimeDerivative
        specific_heat = ${cp2}
        density_name = ${rho2}
        variable = Ts
    []
[]




[AuxKernels]
    [velocity_magnitude]
        type = VectorMagnitudeAux
        variable = velocity_xy
        x = u
        y = v
        # z = w
        execute_on = 'TIMESTEP_BEGIN'
        block = 'coolant'
    []
[]




[Materials]
    [friction_coefficient]
      type = ADGenericFunctorMaterial
      prop_names = 'friction_coefficient'
      prop_values = '0.5'
    []
[]





[FunctorMaterials]
    [functor_constants]
        type = ADGenericFunctorMaterial
        prop_names = 'rho cp k mu drho_dT'
        prop_values = '${rho} ${cp} ${k1} ${mu} ${drho_dT}'
        block = 'coolant'
    []

    [ins_fv]
        type = INSFVEnthalpyFunctorMaterial
        temperature = 'T'
        rho = ${rho}
        block = 'coolant'
    []
    
    [inconel_properties]
        type = ADGenericFunctorMaterial
        prop_names = 'rho2 cp2 k2'
        prop_values = '${rho2} ${cp2} ${k2}'
        block = 'inconel_top inconel_bottom'
    []
[]





[FVBCs]
    [inlet-u]
      type = INSFVInletVelocityBC
      boundary = 'inlet'
      variable = u
      function = ${inlet_velocity}
    []
    [inlet-v]
      type = INSFVInletVelocityBC
      boundary = 'inlet'
      variable = v
      function = 0.001
    []
    # appropriate heat flux condition
    [inlet-T]
      type = FVNeumannBC
      variable = T
      value = '${fparse rho * cp * inlet_velocity * inlet_temp}'
      boundary = 'inlet'
    []
    [no-slip-u]
      type = INSFVNoSlipWallBC
      boundary = 'inner_wall_top inner_wall_bottom'
      variable = u
      function = 0.001
    []
    [no-slip-v]
      type = INSFVNoSlipWallBC
      boundary = 'inner_wall_top inner_wall_bottom'
      variable = v
      function = 0.001
    []
  
    # [symmetry-u]
    #   type = INSFVSymmetryVelocityBC
    #   boundary = 'inner_wall_bottom'
    #   variable = u
    #   u = u
    #   v = v
    #   mu = ${mu}
    #   momentum_component = 'x'
    # []
    # [symmetry-v]
    #   type = INSFVSymmetryVelocityBC
    #   boundary = 'inner_wall_bottom'
    #   variable = v
    #   u = u
    #   v = v
    #   mu = ${mu}
    #   momentum_component = 'y'
    # []
    # [symmetry-p]
    #   type = INSFVSymmetryPressureBC
    #   boundary = 'inner_wall_bottom'
    #   variable = pressure
    # []

    [bottom]
        type = FVDirichletBC
        variable = Ts
        boundary = outer_wall_bottom
        value = 700
    []
    [top]
        type = FVDirichletBC
        variable = Ts
        boundary = outer_wall_top
        value = 700
    []

    [outlet_u]
        type = INSFVMomentumAdvectionOutflowBC
        variable = u
        boundary = 'outlet'
        u = u
        v = v
        momentum_component = 'x'
        rho = ${rho}
    []
    [outlet_v]
        type = INSFVMomentumAdvectionOutflowBC
        variable = v
        boundary = 'outlet'
        u = u
        v = v
        momentum_component = 'y'
        rho = ${rho}
    []

    [outlet_p]
      type = INSFVOutletPressureBC
      boundary = 'outlet'
      variable = pressure
      function = '${outlet_pressure}'
    []


    # [outlet_p]
    #     type = INSFVMassAdvectionOutflowBC
    #     boundary = 'outlet'
    #     variable = pressure
    #     u = u
    #     v = v
    #     rho = ${rho}
    #   []
[]



[FVInterfaceKernels]
    [convection]
        type = FVConvectionCorrelationInterface
        variable1 = T
        variable2 = Ts
        boundary = 'inner_wall_bottom inner_wall_top'
        h = 23241
        T_solid = Ts
        T_fluid = T
        subdomain1 = 'coolant'
        subdomain2 = 'inconel_top inconel_bottom'
        wall_cell_is_bulk = true
      []
[]


[Postprocessors]
    [avg_coolant_temp_pp]
        type = ElementAverageValue
        variable = T
        block = 'coolant'
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [max_coolant_temp_pp]
        type = ElementExtremeValue
        variable = T
        execute_on = 'INITIAL TIMESTEP_END'
        block = 'coolant'
    []
    [avg_channel_temp_pp]
        type = ElementAverageValue
        variable = Ts
        block = 'inconel_bottom'
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [max_channel_temp_pp]
        type = ElementExtremeValue
        variable = Ts
        execute_on = 'INITIAL TIMESTEP_END'
        block = 'inconel_bottom'
    []
    [VEEELOCITY]
        type = ElementExtremeValue
        variable = velocity_xy
        execute_on = 'INITIAL TIMESTEP_END'
        block = 'coolant'
    []
[]


[Executioner]
    type = Transient
    end_time = 300
    dtmin = 1e-8
    dtmax = 0.1
    petsc_options = '-snes_converged_reason -ksp_converged_reason -options_left'
    solve_type = 'NEWTON'
    line_search = 'none'
    nl_max_its = 15
    l_max_its = 500
    nl_rel_tol = 1e-10
    [TimeStepper]
      type = IterationAdaptiveDT
      optimal_iterations = 7
      dt = 0.5
      linear_iteration_ratio = 1e6
      growth_factor = 1.5
    []
[]


[Problem]
    type = NavierStokesProblem
[]


[Outputs]
    file_base = ${output_csv_name}
    exodus = true

    [csv]
        type = CSV
    []

    [console]
        type = Console
        max_rows = 1
    []

    [debug]
        type = MaterialPropertyDebugOutput
    []

    print_linear_residuals = false
    print_nonlinear_residuals = false
    progress = true
[]

[Debug]
    show_material_props = true
[]