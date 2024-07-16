# --- Run Options ---
# Pipe orientation shortcuts
right = '1 0 0'
left = '-1 0 0'
up = '0 1 0'
down = '0 -1 0'

# .csv output name
output_csv_name = 'w1_test'

# Simulation time - "t3" is the FINAL run time
motor_ramp_up_duration = 100.0
motor_ramp_down_duration = 100.0
motor_zero_power = 1000.0
t_end = 30000
t1 = ${motor_ramp_up_duration}
t2 = ${fparse t1 + motor_ramp_down_duration}
t3 = ${fparse t2 + motor_zero_power}
t4 = ${t_end}

# PID options
K_p_value = 7
K_i_value = 0.95
K_d_value = 0.0001

# To be replaced:
temp_in_from_core = 1000
temp_out_at_hx = 300 

# --- Primary Loop ---
# GENERAL
# Initial Conditions
p_m_dot = 9.4                                   # [kg/s]
p_pressure = 9e+6                               # [Pa]
p_initial_vel = 0.00001                         # [m/s]

# MATERIAL PROPERTIES
#   Helium properties - ROUGH
he_gamma = 1.66                                 # = cp/cv
he_molar_mass = 0.004                           # [kg/mol] - Molar mass
he_thermal_conductivity = 0.149                 # [W/m/K]
he_viscosity = 5e-3                             # [Pa*s]

# COMPONENT DETAILS
# Pipe parameters
p_D_pipe = 0.25                                 # [m]
p_A_pipe = ${fparse 0.25 * pi * p_D_pipe^2}     # [m2]

# Pump parameters @Jake-HW-Miles
p_pump_vol_dot = 2                              # [m3/s]
p_pump_head = 350                               # [m]

# --- Secondary Loop ---
# GENERAL
# Initial Conditions
T_cold = 300
T_ambient = 300
p_ambient = 1e5
speed_initial = 0

# COMPONENT DETAILS
# Pipe parameters
s_D_pipe = 0.35
s_A_pipe = ${fparse 0.25 * pi * s_D_pipe^2}

# Motor parameters
I_motor = 1.0
motor_torque_max = 5

# Generator parameters
I_generator = 1.0
#generator_torque_per_shaft_speed = 13.83
generator_torque_per_shaft_speed = 1.383


# Compressor parameters
A_ref_comp = ${s_A_pipe}
V_comp = ${fparse A_ref_comp * 1.0}
I_comp = 1.0

c0_rated_comp = 340.0
rho0_rated_comp = 1.2
c0_t_rated_comp = 670.0
rho0_t_rated_comp = 1.4

eff_comp = 0.79

# Turbine parameters
A_ref_turb = ${s_A_pipe}
V_turb = ${fparse A_ref_turb * 1.0}
I_turb = 1.0

eff_turb = 0.843

# Other (?)
rated_mfr = 20.0

speed_rated_rpm = 9948
speed_rated = ${fparse speed_rated_rpm * 2 * pi / 60.0}

resolution = '${units 10 cm -> m}'

# --- Core calculations ---
mat_reference_temp = 800
n_assemblies = 55

# Coolant rods
cool_rod_per_ass = 18

# Fudged diameters
r_grph = ${fparse 0.040399 * 0.5} # Area based, then made thinner to transfer more heat
r_fuel = 0.01167

# --- Geometry ---
# PRIMARY LOOP (prefix: p_)
p_x1 = 0 
p_y1 = 0
p_pos_1 = '${p_x1} ${p_y1} 0'
p_L1 = 1
p_n_elems1 = ${fparse p_L1 / resolution}

p_x2 = ${fparse p_x1}
p_y2 = ${fparse p_y1 + p_L1}
p_pos_2 = '${p_x2} ${p_y2} 0'
p_L2 = 2
p_n_elems2 = ${fparse p_L2 / resolution}

p_x3 = ${fparse p_x2}
p_y3 = ${fparse p_y2 + p_L2}
p_pos_3 = '${p_x3} ${p_y3} 0'
p_L3 = 1
p_n_elems3 = ${fparse p_L3 / resolution}

p_x4 = ${fparse p_x3}
p_y4 = ${fparse p_y3 + p_L3}
p_pos_4 = '${p_x4} ${p_y4} 0'
p_L4 = 0.5
p_n_elems4 = ${fparse p_L4 / resolution}

p_x4_5 = ${fparse p_x4 + p_L4}
p_y4_5 = ${fparse p_y4}
p_pos_4_5 = '${p_x4_5} ${p_y4_5} 0'
p_L4_5 = 0.5
p_n_elems4_5 = ${fparse p_L4_5 / resolution}

p_x5 = ${fparse p_x4_5 + p_L4_5}
p_y5 = ${fparse p_y4_5}
p_pos_5 = '${p_x5} ${p_y5} 0'
p_L5 = 1
p_n_elems5 = ${fparse p_L5 / resolution}

p_x6 = ${fparse p_x5}
p_y6 = ${fparse p_y5 - p_L5}
p_pos_6 = '${p_x6} ${p_y6} 0'
p_L6 = 3
p_n_elems6 = ${fparse p_L6 / resolution}
# gap_width = 0.002
# hx_wall_thickness = 0.001

# increasing gap width decreases the temperature on secondary side
gap_width = 0.002
hx_wall_thickness = 0.001
p_pos_6_hs = '${fparse p_x6 - gap_width} ${p_y6} 0'

p_x8 = ${fparse p_x6}
p_y8 = ${fparse p_y6 - p_L6}
p_pos_8 = '${p_x8} ${p_y8} 0'
p_L8 = 0.5
p_n_elems8 = ${fparse p_L8 / resolution}

p_x9 = ${fparse p_x8 - p_L8}
p_y9 = ${fparse p_y8}
p_pos_9 = '${p_x9} ${p_y9} 0'
p_L9 = 0.5
p_n_elems9 = ${fparse p_L9 / resolution}

# SECONDARY LOOP (prefix: s_)
s_x1 = ${fparse 1 + gap_width}
s_y1 = 6
s_L1 = 1
s_pos_1 = '${s_x1} ${s_y1} 0'
s_n_elems1 = ${fparse s_L1 / resolution}

s_x2 = ${fparse s_x1 + s_L1}
s_y2 = ${fparse s_y1}
s_L2 = 1
s_pos_2 = '${s_x2} ${s_y2} 0'
s_n_elems2 = ${fparse s_L2 / resolution}

s_x3 = ${fparse s_x2}
s_y3 = ${fparse s_y2 - s_L2}
s_L3 = 5
s_pos_3 = '${s_x3} ${s_y3} 0'
s_n_elems3 = ${fparse s_L3 / resolution}

s_x4 = ${fparse s_x3}
s_y4 = ${fparse s_y3 - s_L3}
s_L4 = 1
s_pos_4 = '${s_x4} ${s_y4} 0'
s_n_elems4 = ${fparse s_L4 / resolution}

s_x5 = ${fparse s_x4 - s_L4}
s_y5 = ${fparse s_y4}
s_L5 = 3
s_pos_5 = '${s_x5} ${s_y5} 0'
s_n_elems5 = ${fparse s_L5 / resolution}

s_x6 = ${fparse s_x5}
s_y6 = ${fparse s_y5 + s_L5}
s_L6 = 2
s_pos_6 = '${s_x6} ${s_y6} 0'
s_n_elems6 = ${fparse s_L6 / resolution}

s_x7 = ${fparse s_x6}
s_y7 = ${fparse s_y6 + s_L6}
s_L7 = 1
s_pos_7 = '${s_x7} ${s_y7} 0'
s_n_elems7 = ${fparse s_L7 / resolution}
s_7+ = ${fparse s_y7 + s_L7 / s_n_elems7}
s_7- = ${fparse s_y7 - s_L6/s_n_elems6}
s_pos_7+ = '${s_x7} ${s_7+} 0'
s_pos_7- = '${s_x7} ${s_7-} 0'

[GlobalParams]
    # Gravity at this stage is constant
    gravity_vector = '0 0 0'

    # Scaling parameters
    rdg_slope_reconstruction = minmod
    # scaling_factor_1phase = '1 1e-2 1e-4'
    # scaling_factor_rhoV = 1
    # scaling_factor_rhouV = 1e-2
    # scaling_factor_rhovV = 1e-2
    # scaling_factor_rhowV = 1e-2
    # scaling_factor_rhoEV = 1e-4

    initial_p = ${p_ambient}
    initial_T = ${T_ambient}
    initial_vel = 0
    initial_vel_x = 0
    initial_vel_y = 0
    initial_vel_z = 0

    fp = fp_air
    closures = thm_closure
    f = 0
[]


[FluidProperties]
    [helium]
        type = IdealGasFluidProperties
        gamma = ${he_gamma}
        k = ${he_thermal_conductivity}
        molar_mass =${he_molar_mass}
        mu = ${he_viscosity}
    []

    [fp_air]
        type = IdealGasFluidProperties
        emit_on_nan = none
    []
[]


[Closures]
    [thm_closure] # Basic closures for 1-phase flow channels
        type = Closures1PhaseTHM
        #    type = Closures1PhaseSimple
    []
[]

[SolidProperties]
    [region1-mat]
      type = ThermalFunctionSolidProperties
      k = 9.8
      cp = 490 
      rho = 8440 
    []
    [fuel-mat]
      type = ThermalFunctionSolidProperties
      k = 16
      cp = 191.67
      rho = 1.4583e4
    []
    [grph-mat]
      type = ThermalFunctionSolidProperties
      k = 24
      cp = 707.7
      rho = 2250
    []
[]

[Components]
    # --- CORE EMULATOR ---
    [core]
        [power]
            type = TotalPower
            power = 15e6                    # Watts = 15MWth
        []

        [reactor]
            type = HeatStructureCylindrical
            position = ${p_pos_2}
            orientation = ${up}
            length = ${p_L2}
            n_elems = 20
            initial_T = 800
            names = 'graph fuel'
            widths = '${r_grph} ${r_fuel}'
            n_part_elems = '2 20'
            solid_properties = 'grph-mat fuel-mat'
            solid_properties_T_ref = '${mat_reference_temp} ${mat_reference_temp}'
            num_rods = ${fparse cool_rod_per_ass * n_assemblies}
        []

        [heat_gen]
            type = HeatSourceFromTotalPower
            hs = core/reactor
            regions = 'fuel'
            power = core/power
            power_shape_function = psf
            power_fraction = 1
        [] 

        [p_core_pipe]
            type = FlowChannel1Phase
            position = ${p_pos_2}
            orientation = ${up}
            length = ${p_L2}
            n_elems = ${p_n_elems2}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_in_from_core}
            initial_vel = ${p_initial_vel}
        []

        [heat_in]
            type = HeatTransferFromHeatStructure1Phase
            flow_channel = core/p_core_pipe
            hs = core/reactor
            hs_side = OUTER
        []
    []

    # --- HEAT EXCHANGER ---
    [hx]
        [p_hx_pipe]
            type = FlowChannel1Phase
            position = ${p_pos_6}
            orientation = ${down}
            length = ${p_L6}
            n_elems = ${p_n_elems6}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [ht1]
            type = HeatTransferFromHeatStructure1Phase
            hs = 'hx/heat_structure'
            hs_side = OUTER
            flow_channel = 'hx/p_hx_pipe'
            Hw = 2e3
            P_hf = ${p_L6}
        []

        [heat_structure]
            type = HeatStructureCylindrical
            position = ${p_pos_6_hs}
            orientation = ${down}
            length = ${p_L6}
            n_elems = ${p_n_elems6}
            names = 'region1'
            n_part_elems = 20
            num_rods = 200000
            solid_properties = 'region1-mat'
            solid_properties_T_ref = '300'
            inner_radius = ${fparse gap_width - hx_wall_thickness}
            widths = ${gap_width}
            initial_T = 300
        []

        [ht2]
            type = HeatTransferFromHeatStructure1Phase
            hs = 'hx/heat_structure'
            hs_side = INNER
            flow_channel = 'hx/s_hx_pipe'
            Hw = 2e3
            P_hf = ${s_L5}
        []

        [s_hx_pipe]
            type = FlowChannel1Phase
            position = '${s_pos_5}'
            length = ${s_L5}
            orientation = ${up}
            n_elems = ${s_n_elems5}
            A = ${s_A_pipe}
        []
    []

    # --- SECONDARY LOOP ---
    [p]
        # Pipes
        [top_pipe_1]
            type = FlowChannel1Phase
            position = ${p_pos_4}
            orientation = ${right}
            length = ${p_L4}
            n_elems = ${p_n_elems4}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [top_pipe_2]
            type = FlowChannel1Phase
            position = ${p_pos_4_5}
            orientation = ${right}
            length = ${p_L4_5}
            n_elems = ${p_n_elems4_5}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [press_pipe]
            type = FlowChannel1Phase
            position = ${p_pos_4_5}
            orientation = ${up}
            length = 0.4
            n_elems = 10
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [down_pipe_1]
            type = FlowChannel1Phase
            position = ${p_pos_5}
            orientation = ${down}
            length = ${p_L5}
            n_elems = ${p_n_elems5}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [bottom_1]
            type = FlowChannel1Phase
            position = ${p_pos_8}
            orientation = ${left}
            length = ${p_L8}
            n_elems = ${p_n_elems8}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [bottom_2]
            type = FlowChannel1Phase
            position = ${p_pos_9}
            orientation = ${left}
            length = ${p_L9}
            n_elems = ${p_n_elems9}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [up_pipe_1]
            type = FlowChannel1Phase
            position = ${p_pos_1}
            orientation = ${up}
            length = ${p_L1}
            n_elems = ${p_n_elems1}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [up_pipe_2]
            type = FlowChannel1Phase
            position = ${p_pos_3}
            orientation = ${up}
            length = ${p_L3}
            n_elems = ${p_n_elems3}
            A = ${p_A_pipe}
            D_h = ${p_D_pipe}
            fp = helium
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel = ${p_initial_vel}
        []
        [bot_pump]
            type = Pump1Phase
            position = ${p_pos_9}
            connections = 'p/bottom_1:out p/bottom_2:in'
            volume = ${p_pump_vol_dot}
            A_ref = ${p_A_pipe}
            head = ${p_pump_head}
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel_x = ${fparse -1*p_initial_vel}
            initial_vel_y = 0
            initial_vel_z = 0
        []
    
        # Junctions
        [j1]
            type = JunctionOneToOne1Phase
            connections = 'p/top_pipe_2:out p/down_pipe_1:in'
        []
        [j2]
            type = JunctionOneToOne1Phase
            connections = 'p/down_pipe_1:out hx/p_hx_pipe:in'
        []
        [j3]
            type = JunctionOneToOne1Phase
            connections = 'hx/p_hx_pipe:out p/bottom_1:in'
        []
        [j5]
            type = JunctionOneToOne1Phase
            connections = 'p/bottom_2:out p/up_pipe_1:in'
        []
        [j6]
            type = JunctionOneToOne1Phase
            connections = 'p/up_pipe_1:out core/p_core_pipe:in'
        []
        [j7]
            type = JunctionOneToOne1Phase
            connections = 'core/p_core_pipe:out p/up_pipe_2:in'
        []
        [j8]
            type = JunctionOneToOne1Phase
            connections = 'p/up_pipe_2:out p/top_pipe_1:in'
        []
        [j9]
            type = VolumeJunction1Phase
            position = ${p_pos_4_5}
            volume = ${fparse 4/3 * pi * (p_D_pipe/2)^3}
            connections = 'p/top_pipe_1:out p/top_pipe_2:in p/press_pipe:in'
            initial_p = ${p_pressure}
            initial_T = ${temp_out_at_hx}
            initial_vel_x = 0
            initial_vel_y = 0
            initial_vel_z = 0
        []
    
        # Pressurizer
        [pressurizer]
            type = InletStagnationPressureTemperature1Phase
            p0 = ${p_pressure}
            T0 = ${temp_out_at_hx}
            input = p/press_pipe:out
        []
    []
    
    [s]
        [pipe1]
            type = FlowChannel1Phase
            position = '${s_pos_1}'
            length = ${s_L1}
            orientation = ${right}
            n_elems = ${s_n_elems1}
            A = ${s_A_pipe}
        []
        [pipe2]
            type = FlowChannel1Phase
            position = '${s_pos_2}'
            length = ${s_L2}
            #length = 1000
            orientation = ${down}
            n_elems = ${s_n_elems2}
            #A = ${s_A_pipe}
            A = 0.00001
        []
        [pipe3]
            type = FlowChannel1Phase
            position = '${s_pos_3}'
            length = ${s_L3}
            orientation = ${down}
            n_elems = ${s_n_elems3}
            A = ${s_A_pipe}
        []
        [pipe4]
            type = FlowChannel1Phase
            position = '${s_pos_4}'
            length = ${s_L4}
            orientation = ${left}
            n_elems = ${s_n_elems4}
            A = ${s_A_pipe}
        []
        [pipe6]
            type = FlowChannel1Phase
            position = '${s_pos_6}'
            length = ${s_L6}
            orientation = ${up}
            n_elems = ${s_n_elems6}
            A = ${s_A_pipe}
        []
        [pipe7]
            type = FlowChannel1Phase
            position = '${s_pos_7}'
            length = ${s_L7}
            orientation = ${up}
            n_elems = ${fparse s_n_elems7}
            A = ${s_A_pipe}
        []
        # Junctions
        [junction1]
            type = JunctionOneToOne1Phase
            connections = 's/pipe1:out s/pipe2:in'
        []
        [junction2]
            type = JunctionOneToOne1Phase
            connections = 's/pipe3:out s/pipe4:in'
        []
        [junction3]
            type = JunctionOneToOne1Phase
            connections = 's/pipe4:out hx/s_hx_pipe:in'
        []
        [junction4]
            type = JunctionOneToOne1Phase
            connections = 'hx/s_hx_pipe:out s/pipe6:in'
        []
        [junction5]
            type = JunctionOneToOne1Phase
            connections = 's/pipe7:out s/pipe1:in'
        []

        [ultimate_heat_sink]
            type = HeatTransferFromSpecifiedTemperature1Phase
            flow_channel = s/pipe2
            T_wall = ${T_cold}
            # Hw = 100000000
        []

        # [heat_structure]
        #     type = HeatStructurePlate
        #     position = ${s_pos_3}
        #     orientation = ${down}
        #     length = ${s_L3}
        #     n_elems = ${s_n_elems3}
        #     depth = 1
        #     names = 'region1'
        #     solid_properties = 'region1-mat'
        #     solid_properties_T_ref = '300'
        #     widths = '1'
        #     n_part_elems = '10'
        #     initial_T = ${T_cold}
        #   []
        # [temp_outside]
        #     type = HSBoundarySpecifiedTemperature
        #     hs = s/heat_structure
        #     boundary = s/heat_structure:out
        #     T = ${T_cold}
        # []
    []

    # [ht1]
    #     type = HeatTransferFromHeatStructure1Phase
    #     hs = 'heat_structure'
    #     hs_side = INNER
    #     flow_channel = 's/pipe2'
    #     Hw = 10000
    #     P_hf = ${p_L6}
    # []
    # [heat_structure]
    #     type = HeatStructureCylindrical
    #     position = ${s_pos_2}
    #     orientation = ${down}
    #     length = ${s_L2}
    #     n_elems = ${s_n_elems2}
    #     names = 'region1'
    #     n_part_elems = ${s_n_elems2}
    #     num_rods = 200
    #     solid_properties = 'region1-mat'
    #     solid_properties_T_ref = '300'
    #     inner_radius = 0.001
    #     widths = 0.01
    #     initial_T = 300
    # []
    # [ambient_convection]
    #     type = HSBoundaryAmbientConvection
    #     boundary = 'heat_structure:outer'
    #     hs = heat_structure
    #     T_ambient = ${T_cold}
    #     htc_ambient = 1e10
    # []

    [shaft]
        type = Shaft
        connected_components = 'motor compressor turbine generator'
        initial_speed = ${speed_initial}
    [] 
    [motor]
        type = ShaftConnectedMotor
        inertia = ${I_motor}
        torque = 0 # controlled
    []
    [generator]
        type = ShaftConnectedMotor
        inertia = ${I_generator}
        torque = generator_torque_fn
    []
    [compressor]
        type = ShaftConnectedCompressor1Phase
        position = '${s_pos_3}'
        inlet = 's/pipe2:out'
        outlet = 's/pipe3:in'
        A_ref = ${A_ref_comp}
        volume = ${V_comp}
        omega_rated = ${speed_rated}
        mdot_rated = ${rated_mfr}
        c0_rated = ${c0_rated_comp}
        rho0_rated = ${rho0_rated_comp}
        speeds = '0.5208 0.6250 0.7292 0.8333 0.9375'
        Rp_functions = 'rp_comp1 rp_comp2 rp_comp3 rp_comp4 rp_comp5'
        eff_functions = 'eff_comp1 eff_comp2 eff_comp3 eff_comp4 eff_comp5'
        min_pressure_ratio = 1
        speed_cr_I = 0
        inertia_const = ${I_comp}
        inertia_coeff = '${I_comp} 0 0 0'
        speed_cr_fr = 0
        tau_fr_const = 0
        tau_fr_coeff = '0 0 0 0'
    []
    [turbine]
        type = ShaftConnectedCompressor1Phase
        position = '${s_pos_7}'
        inlet = 's/pipe6:out'
        outlet = 's/pipe7:in'
        A_ref = ${A_ref_turb}
        volume = ${V_turb}
        treat_as_turbine = true
        omega_rated = ${speed_rated}
        mdot_rated = ${rated_mfr}
        c0_rated = ${c0_t_rated_comp}
        rho0_rated = ${rho0_t_rated_comp}
        speeds = '0 0.5208 0.6250 0.7292 0.8333 0.9375'
        Rp_functions = 'rp_turb0 rp_turb1 rp_turb2 rp_turb3 rp_turb4 rp_turb5'
        eff_functions = 'eff_turb1 eff_turb1 eff_turb2 eff_turb3 eff_turb4 eff_turb5'
        min_pressure_ratio = 1
        speed_cr_I = 0
        inertia_const = ${I_turb}
        inertia_coeff = '${I_turb} 0 0 0'
        speed_cr_fr = 0
        tau_fr_const = 0
        tau_fr_coeff = '0 0 0 0'
    []
[]


[ControlLogic]
# Three consectutive control blocks used to set the control of the secondary loop pump:
    [mass_flow_rate_set_point] # Defines the TARGET VALUE for mass flow rate in pump
        type = GetFunctionValueControl
        function = ${p_m_dot}
    []

    [pid_for_secondary_pump]
        type = PIDControl
        initial_value = 0
        set_point = mass_flow_rate_set_point:value
        input = mass_flow_rate_in_pump
        K_p = ${K_p_value}
        K_i = ${K_i_value}
        K_d = ${K_d_value}
    []

    [set_p_pump_head]
        type = SetComponentRealValueControl
        component = p/bot_pump
        parameter = head
        value = pid_for_secondary_pump:output
    []

    [motor_ctrl]
        type = TimeFunctionComponentControl
        component = motor
        parameter = torque
        function = motor_torque_fn
      []
[]

[Postprocessors]
    
    [p_ratio_turb]
        type = ParsedPostprocessor
        pp_names = 'p_in_turb p_out_turb'
        expression = 'p_in_turb / p_out_turb'
        execute_on = 'INITIAL TIMESTEP_END'
      []
    
    [set_point_tracker]
        type = RealControlDataValuePostprocessor
        control_data_name = mass_flow_rate_set_point:value
    []

    [mass_flow_rate_in_pump]
        type = ADFlowJunctionFlux1Phase
        junction = p/j3
        connection_index = 1
        boundary = p/bottom_1:in
        equation = mass
    []

    [p_in_turb]
        type = PointValue
        variable = p
        point = ${s_pos_7-}
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [p_out_turb]
        type = PointValue
        variable = p
        point = ${s_pos_7+}
        execute_on = 'INITIAL TIMESTEP_END'
    []
    
    [T_in_hx]
        type = PointValue
        variable = T
        point = ${s_pos_5}
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [T_out_hx]
        type = PointValue
        variable = T
        point = ${s_pos_6}
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [T_out_HS]
        type = PointValue
        variable = T
        point = ${s_pos_3}
        execute_on = 'INITIAL TIMESTEP_END'
    []

    [T_in_core]
        type = PointValue
        variable = T
        point = ${p_pos_3}
        execute_on = 'INITIAL TIMESTEP_END'
    []

    [T_out_core]
        type = PointValue
        variable = T
        point = ${p_pos_4}
        execute_on = 'INITIAL TIMESTEP_END'
    []

    [v_]
        type = PointValue
        variable = T
        point = ${s_pos_3}
        execute_on = 'INITIAL TIMESTEP_END'
    []

    [error_tracker]
        type = DifferencePostprocessor
        value1 = set_point_tracker
        value2 = mass_flow_rate_in_pump
    []

    [heating_rate]
        type = ADHeatRateConvection1Phase
        block = 'hx/s_hx_pipe'
        T = T
        T_wall = T_wall
        Hw = Hw
        P_hf = P_hf
        execute_on = 'INITIAL TIMESTEP_END'
    []
    # [cooling_rate]
    #     type = ADHeatRateConvection1Phase
    #     block = 's/pipe3'
    #     T = T
    #     T_wall = T_wall
    #     Hw = Hw
    #     P_hf = P_hf
    #     execute_on = 'INITIAL TIMESTEP_END'
    # []
    [motor_torque]
        type = RealComponentParameterValuePostprocessor
        component = motor
        parameter = torque
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [motor_power]
        type = FunctionValuePostprocessor
        function = motor_power_fn
        execute_on = 'INITIAL TIMESTEP_END'
        indirect_dependencies = 'motor_torque shaft:omega'
    []
    [generator_torque]
        type = ShaftConnectedComponentPostprocessor
        quantity = torque
        shaft_connected_component_uo = generator:shaftconnected_uo
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [generator_power]
        type = FunctionValuePostprocessor
        function = generator_power_fn
        execute_on = 'INITIAL TIMESTEP_END'
        indirect_dependencies = 'generator_torque shaft:omega'
    []
    [shaft_speed]
        type = ScalarVariable
        variable = 'shaft:omega'
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [mfr_comp]
        type = ADFlowJunctionFlux1Phase
        boundary = s/pipe2:out
        connection_index = 0
        equation = mass
        junction = compressor
    []
    [mfr_turb]
        type = ADFlowJunctionFlux1Phase
        boundary = s/pipe6:out
        connection_index = 0
        equation = mass
        junction = turbine
    []
[]

[Preconditioning]
    [pc]
        type = SMP
        full = true
    []
[]
  
[Executioner]
    type = Transient
    start_time = 0
    scheme = 'bdf2'

    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 0.001
    []

    dtmax = 2
    end_time = ${t4}
  
    line_search = basic
    solve_type = NEWTON
  
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
  
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-8
    nl_max_its = 20

    l_tol = 1e-4
    l_max_its = 40  
[]
  
[Outputs]
    file_base = ${output_csv_name}

    exodus = true
    print_linear_residuals = false
    print_nonlinear_residuals = false
    
    [csv]
        type = CSV
    []

    [console]
        type = Console
        max_rows = 5
        show = 'generator_power shaft_speed p_ratio_turb T_in_hx T_out_hx T_out_HS T_in_core T_out_core'
    []
[]

[Functions]
    [psf] # Power Shape Function (PSF) -> Used for core
      type = ParsedFunction
      expression = 1
    []
    [motor_torque_fn]
      type = PiecewiseLinear
      x = '0 ${t1} ${t2} ${t3}'
      y = '0 ${motor_torque_max} 1 1'
    []
    [motor_power_fn]
      type = ParsedFunction
      expression = 'torque * speed'
      symbol_names = 'torque speed'
      symbol_values = 'motor_torque shaft:omega'
    []
    [generator_torque_fn]
      type = ParsedFunction
      expression = 'slope * t'
      symbol_names = 'slope'
      symbol_values = '${generator_torque_per_shaft_speed}'
    []
    [generator_power_fn]
      type = ParsedFunction
      expression = 'torque * speed'
      symbol_names = 'torque speed'
      symbol_values = 'generator_torque shaft:omega'
    []
    [htc_wall_fn]
      type = PiecewiseLinear
        x = '0 ${t1} ${t2} ${t_end}'
        y = '0 1e3 1e3 1e3'
    []

    [rp_comp1]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_comp1.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_comp2]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_comp2.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_comp3]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_comp3.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_comp4]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_comp4.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_comp5]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_comp5.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_turb0]
        type = ConstantFunction
        value = 1
      []
      [rp_turb1]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_turb1.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_turb2]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_turb2.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_turb3]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_turb3.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_turb4]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_turb4.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
      [rp_turb5]
        type = PiecewiseLinear
        data_file = 'csv_store_20MFR/rp_turb5.csv'
        x_index_in_file = 0
        y_index_in_file = 1
        format = columns
        extrap = true
      []
    
      # compressor efficiency
      [eff_comp1]
        type = ConstantFunction
        value = ${eff_comp}
      []
      [eff_comp2]
        type = ConstantFunction
        value = ${eff_comp}
      []
      [eff_comp3]
        type = ConstantFunction
        value = ${eff_comp}
      []
      [eff_comp4]
        type = ConstantFunction
        value = ${eff_comp}
      []
      [eff_comp5]
        type = ConstantFunction
        value = ${eff_comp}
      []
    
      # turbine efficiency
      [eff_turb1]
        type = ConstantFunction
        value = ${eff_turb}
      []
      [eff_turb2]
        type = ConstantFunction
        value = ${eff_turb}
      []
      [eff_turb3]
        type = ConstantFunction
        value = ${eff_turb}
      []
      [eff_turb4]
        type = ConstantFunction
        value = ${eff_turb}
      []
      [eff_turb5]
        type = ConstantFunction
        value = ${eff_turb}
      []
[]