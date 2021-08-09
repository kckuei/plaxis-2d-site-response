import os
#import pickle
import numpy as np
import pandas as pd
from shapely.geometry import (MultiLineString, LineString, MultiPolygon, 
                              Polygon, GeometryCollection, Point)
# from shapely.ops import split
import matplotlib.pylab as plt
from matplotlib import cm
from scipy import interpolate

import utils as ut


#-----------------------------------------------------------------------------
# Meshing
#-----------------------------------------------------------------------------
# Load polyzone and surface geometries
UserGeometry, _ = ut.read_excel_poly("User Inputs/User_Geometry.xlsx") 


GenerateZones = ['Sand', 'Clay', 'Till']  # Specify zones to generate

NxNy = {#'Sand':(40, 20),                  # Specify number of columns and rows
        'Sand':(1, 1),
        'Clay':(1, 1), 
        'Till':(1, 1),
        }

# Create the polyzone mesh
PolyMeshCollections, PolyCoordCollections = ut.mesh_grid(UserGeometry,
                                                        GenerateZones,
                                                        NxNy,
                                                        case="topo",
                                                        plot_mesh=False)

#-----------------------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------------------
# Define soil parameter mapping functions for zones of interest

def sand_props():  
    #-------------------------------------------------------------------------
    # User-defined mapping
    #-------------------------------------------------------------------------
    zone = 'Sand'
    model = 'PM4Sand'
    
    # DR (random field)
    realiz_idx = 0
    start_elev = 8
    load_RF_path = r"User Inputs\Random Field\17-02-2021_06-01AM\Realizations_EW.pkl"
    DR, __, __ = ut.map_DR(load_RF_path, realiz_idx, start_elev,
                                   PolyMeshCollections, 
                                   PolyCoordCollections, 
                                   zone, 
                                   DR_min = 0.35, DR_max = 0.85,
                                   preview=False)
    # Modulus
    Go = 61.1 + 1088.1*DR 
    
    # Contraction parameter
    df = pd.read_csv(r"User Inputs\User_DR-hpo Calibration_IB08.csv")
    f_hpo = interpolate.interp1d(df.DR, df.hpo)
    hpo = f_hpo(DR)
    
    # Unit Weight
    unit_weight = 20
    
    # Permeability
    perm = 1e-4

    # For effective stress calculations
    unit_weight = 18
    unit_weight_water = 9.81
    
    # Ko, and poisson ratio
    Ko = 0.5
    nu = min(0.45, Ko/(1+Ko))
    
    # Interpolating functions
    # ..Surface
    f_surf = interpolate.interp1d(UserGeometry['Surface'].x,
                                  UserGeometry['Surface'].y)
    # ..Water level
    f_water = interpolate.interp1d(UserGeometry['Water Level'].x,
                                   UserGeometry['Water Level'].y)
    
    # Effective stress and elastic moduli
    N = len(PolyCoordCollections[zone])
    xc, yc = ut.get_zone_xyc(PolyMeshCollections[zone])
    Gref     = np.zeros(N)
    Eref     = np.zeros(N)
    for i, poly in enumerate(PolyMeshCollections[zone]):
    
        # Effective stresses
        water_elev = f_water(xc[i])        
        ysurf = f_surf(xc[i])
        sv = (ysurf - yc[i])*unit_weight + max(0, (water_elev-ysurf)*unit_weight_water)
        uo = max(0, (water_elev - yc[i])*unit_weight_water)
        esv = sv - uo
        
        # Elastic moduli
        p = esv*(1+Ko)/3
        Gref[i] = Go[i]*101.3*(p/101.3)**0.5
        Eref[i] = 2*Gref[i]*(1+nu)
        
        # Report any negative stresses
        if p < 0:
            print("NEGATIVE EFFECTIVE STRESSES COMPUTED")
            print(f"yc={yc[i]}, ysurf={ysurf}, waterelev={water_elev}")
            print(f"Gref:{Gref[i]}, sv:{sv}, uo:{uo}, esv:{esv}, p:{p}")
    

    
    #-------------------------------------------------------------------------
    # Pack parameters into dictionary
    #-------------------------------------------------------------------------
    SoilParamCollections = []
    N = len(PolyCoordCollections[zone])
    for i in range(N): 
        
        params_dict = {'zone':  zone,
                       'model': model,
                       'xc':    xc[i],
                       'yc':    yc[i],
                       'nu':    round(nu, 2),
                       'Gref':  round(Gref[i], 2),
                       'Eref':  round(Eref[i], 2),
                       'DR':    round(DR[i], 2),
                       'Go':    round(Go[i], 2),
                       'hpo':   round(hpo[i] ,2),
                       'k':     perm,
                       'unit_weight': unit_weight,
                       'rayleigh_alpha': 0,          # Temporary
                       'rayleigh_beta': 0            # Temporary
                      }
        SoilParamCollections += [params_dict]            
    return SoilParamCollections

def clay_props():
    #-------------------------------------------------------------------------
    # User-defined mapping
    #-------------------------------------------------------------------------
    zone = 'Clay'
    model = 'PM4Silt'
    
    # For effective stress calculations
    unit_weight = 18
    unit_weight_water = 9.81
    water_elev = 2.5
    
    # Permeability
    perm = 1e-8
    
    # Interpolating functions
    # ..Top of clay
    f_ytopclay = interpolate.interp1d(UserGeometry['Clay Top'].x,
                                      UserGeometry['Clay Top'].y)
    # ..Surface
    f_surf = interpolate.interp1d(UserGeometry['Surface'].x,
                                  UserGeometry['Surface'].y)
    # ..Water level
    f_water = interpolate.interp1d(UserGeometry['Water Level'].x,
                                   UserGeometry['Water Level'].y)
    
    # Main calibration parameters (Su_ratio, Go, hpo)
    N = len(PolyMeshCollections[zone])
    xc, yc = ut.get_zone_xyc(PolyMeshCollections[zone])
    Su_ratio = np.zeros(N)
    Go       = np.zeros(N)
    hpo      = np.zeros(N)
    Gref     = np.zeros(N)
    Eref     = np.zeros(N)
    for i, poly in enumerate(PolyMeshCollections[zone]):
    
        ytop_clay = f_ytopclay(xc[i])
        zc = ytop_clay - yc[i]
        
        # Stress history, Ko and poisson ratio
        OCR = max(1, 1.15 - zc/1000)
        Ko = 0.5*OCR
        nu = min(0.45, Ko/(1+Ko))

        # Su ratio
        Su_ratio[i] = 0.22*OCR**0.8
        Su_ratio[i] = max(0.25, min(1.0,Su_ratio[i])) 
    
        # Effective stresses
        ysurf = f_surf(xc[i])
        water_elev = f_water(xc[i])
        sv = (ysurf - yc[i])*unit_weight + max(0, (water_elev-ysurf)*unit_weight_water)
        uo = max(0, (water_elev - yc[i])*unit_weight_water)
        esv = sv - uo
    
        # Elastic and NL moduli       
        Su = Su_ratio[i]*esv
        Gref[i] = 915*(OCR**(-0.343))*Su
        Go[i] = Gref[i]/(101.3*(esv/101.3)**0.5)
        Eref[i] = 2*Gref[i]*(1+nu)

        
        # Contraction parameter
        if Su_ratio[i] >= 0.5:        
            hpo[i] = 65.473 + 21.839*np.log(Su_ratio[i])
        else:
            hpo[i] = 80.0+43.281*np.log(Su_ratio[i])
            
        # Report any negative stresses
        if esv < 0:
            print("NEGATIVE EFFECTIVE STRESSES COMPUTED")
            print(f"yc={yc[i]}")
            print(f"Su_ratio:{Su_ratio}, sv:{sv}, uo:{uo}, esv:{esv}")
        
        
    #-------------------------------------------------------------------------
    # Pack parameters into dictionary
    #-------------------------------------------------------------------------
    SoilParamCollections = []
    N = len(PolyCoordCollections[zone])
    for i in range(N): 
        
        params_dict = {'zone':      zone,
                       'model':     model,
                       'xc':        xc[i],
                       'yc':        yc[i],
                       'nu':        round(nu, 2), 
                       'Gref':      round(Gref[i], 2),
                       'Eref':      round(Eref[i], 2),
                       'Su_ratio':  round(Su_ratio[i],2),
                       'Go':        round(Go[i],2),
                       'hpo':       round(hpo[i],2),
                       'k':         perm,
                       'unit_weight': unit_weight,
                       'rayleigh_alpha': 0,          # Temporary
                       'rayleigh_beta': 0            # Temporary
                      }
        SoilParamCollections += [params_dict]            
    return SoilParamCollections

def till_props():
    #-------------------------------------------------------------------------
    # User-defined mapping
    #-------------------------------------------------------------------------
    zone = 'Till'
    model = 'Elastic'
    
    # Unit weight
    unit_weight = 21
    
    # Permeability
    perm = 1e-12
    
    # Elastic parameters
    Vs = 460 
    nu = 0.25
    Gref = (unit_weight/9.81)*Vs**2
    Eref = 2*Gref*(1+nu)
    
    xc, yc = ut.get_zone_xyc(PolyMeshCollections[zone])
            
    #-------------------------------------------------------------------------
    # Pack parameters into dictionary
    #-------------------------------------------------------------------------
    SoilParamCollections = []
    N = len(PolyCoordCollections[zone])
    for i in range(N): 
        
        params_dict = {'zone':      zone,
                       'model':     model,
                       'xc':        xc[i],
                       'yc':        yc[i],
                       'Vs':        round(Vs, 2),
                       'nu':        round(nu, 2), 
                       'Gref':      round(Gref, 2),
                       'Eref':      round(Eref, 2),
                       'k':         perm,
                       'unit_weight': unit_weight,
                       'rayleigh_alpha': 0,          # Temporary
                       'rayleigh_beta': 0            # Temporary
                      }
        SoilParamCollections += [params_dict]            
    return SoilParamCollections
    return SoilParamCollections

# Create soil parameters collection
SoilParamCollections = {'Sand': sand_props(),
                        'Clay': clay_props(),
                        'Till': till_props()}

#-----------------------------------------------------------------------------
# PLAXIS Testing
#-----------------------------------------------------------------------------

# Model properties
ModelSetup = {"Title" : "",
              "Company" : "Golder Associates Ltd.",
              "UnitLength" : "m",
              "UnitForce" : "kN",
              "UnitTime" : "day",
              "ElementType": "6-Noded",
              "ModelType": "Plane strain",
              "WaterWeight" : 10.0,
              "Comments" : ""
             }

# Model boundaries
ModelBoundaries = {'xmin' : 0,  
                   'xmax' : 2610,
                   'ymin' : -120,
                   'ymax' :  20,
                  }


# Base BC and 
BaseBC = 'compliant'    # 'compliant' for outcrop, 'fixed' for within

# Motions
Motion = []
motion_path = r'User Inputs/Motions A2475/DEEPSOIL Acc'
for filename in os.listdir(motion_path):
    if filename.endswith('.txt'):
        fname = os.path.join(motion_path, filename)
        Motion += [pd.read_csv(fname , skiprows = 1, delimiter = '\t', 
                               names = ['Time (s)','Accel (g)'])]
MotionNames = [x.split('_')[0] for x in os.listdir(motion_path)]

#-----------------------------------------------------------------------------
# PLAXIS connection
#-----------------------------------------------------------------------------
from plxscripting.easy import *
s_i, g_i = new_server('localhost', 10000, password='s=5u+SkvSwL4u>zD') 
s_o, g_o = new_server('localhost', 10001, password='s=5u+SkvSwL4u>zD')  

s_i.new()

#-----------------------------------------------------------------------------
# PLAXIS Model Setup
#-----------------------------------------------------------------------------

# Define project properties
g_i.setproperties(*tuple(zip(list(ModelSetup.keys()),
                             list(ModelSetup.values()))))
        
# Sets up model domain
g_i.SoilContour.initializerectangular(ModelBoundaries['xmin'],
                                      ModelBoundaries['ymin'],
                                      ModelBoundaries['xmax'],
                                      ModelBoundaries['ymax'])

#-----------------------------------------------------------------------------
# PLAXIS Geometry/ Polyzones
#-----------------------------------------------------------------------------
g_i.gotostructures()

plx_polys = []
plx_soils = []
for zone in PolyMeshCollections.keys():
    
    Nzones = len(PolyMeshCollections[zone])
    for i in range(Nzones):
        
        points = PolyMeshCollections[zone][i].exterior.coords
        Poly, Soil = g_i.polygon(*points)
        
        plx_polys += [Poly]
        plx_soils += [Soil]
    
#-----------------------------------------------------------------------------
# PLAXIS Material Creation
#-----------------------------------------------------------------------------
# DR color binning
DR_binning_dict = {'DR_min' : 0.25,
                   'DR_max' : 0.85,
                   'DR_increment' : 0.05,
                  }

DR_bins = np.arange(DR_binning_dict['DR_min'], 
                    DR_binning_dict['DR_max'] + 
                    DR_binning_dict['DR_increment'], 
                    DR_binning_dict['DR_increment'])
cmap = [cm.Reds(x) for x in np.linspace(0,1.0,len(DR_bins))]

DR_colors = []
for i, DR in enumerate(DR_bins):
    
    R = int(cmap[i][0]*255)
    G = int(cmap[i][1]*255)
    B = int(cmap[i][2]*255)
    color = ut.get_RGB_number(R, G, B)
    DR_colors.append(color)


def create_mat_elastic(g_i, Phase, soil_info):
    
    name = "Soil_" + str(soil_info['zone']) + '_Elastic_'+ Phase
    drainage_type = 0 #drained
    
    color = 8421504
   
    soil_params = [("MaterialName", name),
                    ("SoilModel", 1),
                    ("Colour", color),
                    ("DrainageType", drainage_type),
                    ("gammaUnsat", soil_info['unit_weight']),
                    ("gammaSat", soil_info['unit_weight']),
                    ("perm_primary_horizontal_axis", 4000),
                    ("perm_vertical_axis", 4000),
                    ("Eref", soil_info['Eref']),
                    ("nu", soil_info['nu']),
                    ("RayleighAlpha", soil_info['rayleigh_alpha']),
                    ("RayleighBeta", soil_info['rayleigh_beta']),
                    ]
    
    return g_i.soilmat(*soil_params)

def create_mat_pm4sand(g_i, Phase, soil_info, end_of_shaking=False):
    
    name = "Soil_" + str(soil_info['zone']) + '_PM4Sand_' + Phase
    drainage_type = 1 #undrainedA - Dynamic with Conlsolidation calculation type only
    
    
    idx = np.abs(DR_bins - soil_info['DR']).argmin()
    color = DR_colors[idx]
    
    ## Add min/max DR bounds check later
    
    if end_of_shaking:
        kh = 4000
        kv = 4000
        post_shake = 1
    else:
        mpsec_to_mpday = 86400
        kh = soil_info['k']*mpsec_to_mpday
        kv = soil_info['k']*mpsec_to_mpday
        post_shake = 0
    
    soil_params = [("MaterialName", name),
                    ("SoilModel", 100),
                    ("UserDLLName", "pm4sand64.dll"),
                    ("UserModel", "PM4Sand"),
                    ("Colour", color),
                    ("DrainageType", drainage_type),
                    ("gammaUnsat", soil_info['unit_weight']),
                    ("gammaSat", soil_info['unit_weight']),
                    ("perm_primary_horizontal_axis", kh),
                    ("perm_vertical_axis", kv),
                    ("PsiUnsat", 0.5),
                    ("User1", soil_info['DR']),     #Dr
                    ("User2", soil_info['Go']),     #G0
                    ("User3", soil_info['hpo']),    #hp0
                    ("User4", 101.3),               #pA
                    ("User5", 0.8),                 #e_max
                    ("User6", 0.5),                 #e_min
                    ("User7", 0.5),                 #nb
                    ("User8", 0.1),                 #nd
                    ("User9", 33),                  #phi_cv
                    ("User10", 0.3),                #nu
                    ("User11", 10),                 #Q
                    ("User12", 1.5),                #R
                    ("User13", post_shake),         #PostShake
                    ("RayleighAlpha", soil_info['rayleigh_alpha']),
                    ("RayleighBeta", soil_info['rayleigh_beta']),                      
                    ("phi", 2),
                    ("Gref", 1000),         
                    ]
    
    return g_i.soilmat(*soil_params)

def create_mat_pm4silt(g_i, Phase, soil_info, end_of_shaking=False):
    
    name = "Soil_" + str(soil_info['zone']) + '_PM4Silt_' + Phase
    drainage_type = 1 #undrainedA - Dynamic with Conlsolidation calculation type only
    
    color = 11750174
    
    if end_of_shaking:
        kh = 4000
        kv = 4000
        post_shake = 1
    else:
        mpsec_to_mpday = 86400
        kh = soil_info['k']*mpsec_to_mpday
        kv = soil_info['k']*mpsec_to_mpday
        post_shake = 0
    
    soil_params = [("MaterialName", name),
                    ("SoilModel", 100),
                    ("UserDLLName", "pm4silt64.dll"),
                    ("UserModel", "PM4Silt"),
                    ("Colour", color),
                    ("DrainageType", drainage_type),
                    ("gammaUnsat", soil_info['unit_weight']),
                    ("gammaSat", soil_info['unit_weight']),
                    ("perm_primary_horizontal_axis", kh),
                    ("perm_vertical_axis", kv),
                    ("PsiUnsat", 0.5),
                    ("User1", soil_info['Su_ratio']),           #Su_ratio
                    #("User2", soil_info['Su']),                #Su
                    ("User3", soil_info['Go']),                 #Go
                    ("User4", soil_info['hpo']),                #hpo
                    ("User5", 101.3),           #Patm
                    ("User6", 0.75),            #nG
                    ("User7", 0.5),             #ho
                    ("User8", 0.90),            #eo
                    ("User9", 0.060),           #lambda
                    ("User10", 32.0),           #phi_cv    ?? Can't assign ??
                    ("User11", 0.80),           #nb,wet
                    ("User12", 0.5),            #nb,dry
                    ("User13", 0.30),           #nd
                    ("User14", 0.80),           #Ado
                    ("User15", 0),              #ru,max     fcn, 0 for default
                    ("User16", 0),              #zmax       fcn, 0 for default
                    ("User17", 100),            #Cz
                    ("User18", 0),              #Cepsilon   fcn, 0 for default
                    ("User19", 3),              #CGD
                    ("User20", 4),              #Ckalphaf
                    ("User21", 0.30),           #nu
                    ("User22", post_shake),     #PostShake
                    ("User23", 2.0),            #CGconsol
                    ("User24", 1.0),            #FSu
                    ("User25", 0),              #psi_R      Alternative.. 
                    ("User26", 0),              #dR         ..two-stage params
                    ("User27", 0),              #phi_c      ..for variable Su
                    ("User28", 0),              #Cc         ..
                    ("RayleighAlpha", soil_info['rayleigh_alpha']),
                    ("RayleighBeta", soil_info['rayleigh_beta']),                      
                    ("phi", 2),
                    ("Gref", 1000),         
                    ]
    return g_i.soilmat(*soil_params)


Mat_Elastic = []
Mat_Dynamic = []
Mat_Postliq = []
for zone in SoilParamCollections.keys():
    
    for i, soil_info in enumerate(SoilParamCollections[zone]):
        
        if soil_info['model'] == 'Elastic':
            MatSet_Elastic = create_mat_elastic(g_i, f'{i}_Elastic', soil_info)
            MatSet_Dynamic = create_mat_elastic(g_i, f'{i}_Dynamic', soil_info)
            MatSet_Postliq = create_mat_elastic(g_i, f'{i}_PostLiq', soil_info)
            
        if soil_info['model'] == 'PM4Sand':
            MatSet_Elastic = create_mat_elastic(g_i, f'{i}_Elastic', soil_info)
            MatSet_Dynamic = create_mat_pm4sand(g_i, f'{i}_Dynamic', soil_info)
            MatSet_Postliq = create_mat_pm4sand(g_i, f'{i}_PostLiq', soil_info,
                                                      end_of_shaking=True)
            
        if soil_info['model'] == 'PM4Silt':
            MatSet_Elastic = create_mat_elastic(g_i, f'{i}_Elastic', soil_info)
            MatSet_Dynamic = create_mat_pm4silt(g_i, f'{i}_Dynamic', soil_info)
            MatSet_Postliq = create_mat_pm4silt(g_i, f'{i}_PostLiq', soil_info,
                                                      end_of_shaking=True)
    
        Mat_Elastic.append(MatSet_Elastic)
        Mat_Dynamic.append(MatSet_Dynamic)
        Mat_Postliq.append(MatSet_Postliq)


# Set elastic material for initial stage
for i in range(len(Mat_Elastic)):
     g_i.Soils[i].Material = Mat_Elastic[i]

#-----------------------------------------------------------------------------
# PLAXIS Stagd Construction / Phases
#-----------------------------------------------------------------------------
g_i.gotostages()

# Initialize
N_motion = len(Motion)
phase_dynamic = [None]*N_motion
phase_postliq = [None]*N_motion

# Create head
phase_initial_stress = g_i.phase(g_i.InitialPhase)
g_i.activate(g_i.Polygons,g_i.InitialPhase)
N_phase_initial = len(g_i.Phases)

# Create dynamic branches
for count in range(N_motion):
    phase_dynamic[count] = g_i.phase(phase_initial_stress)
    phase_dynamic[count].Identification = MotionNames[count] + '_Dynamic'

# Create postliq branches
for count in range(N_motion):
    phase_postliq[count] = g_i.phase(phase_dynamic[count])
    phase_postliq[count].Identification = MotionNames[count] + '_PostLiq'

#[phase.Identification.value for phase in g_i.Phases]

# Assign dynamic material sets
dyn_branch_start = N_phase_initial
for i in range(len(g_i.Soils)):    
    for phase_dyn in g_i.Phases[dyn_branch_start:dyn_branch_start + N_motion]:
        g_i.Soils[i].Material[phase_dyn] = Mat_Dynamic[i]

# Assign postliq material sets
post_branch_start = dyn_branch_start + N_motion
for i in range(len(g_i.Soils)):    
    for phase_post in g_i.Phases[post_branch_start:post_branch_start + N_motion]:
        g_i.Soils[i].Material[phase_post] = Mat_Postliq[i]

g_i.gotostages()
g_i.activate(g_i.Polygons, g_i.InitialPhase)
#-----------------------------------------------------------------------------
# PLAXIS Water Level
#-----------------------------------------------------------------------------
g_i.gotoflow()
g_i.deactivate(g_i.Water, g_i.InitialPhase)


# Interpolating functions
# ..Surface
f_surf = interpolate.interp1d(UserGeometry['Surface'].x,
                              UserGeometry['Surface'].y)
# ..Water level
f_water = interpolate.interp1d(UserGeometry['Water Level'].x,
                               UserGeometry['Water Level'].y)


global_water_level = g_i.waterlevel(*zip(UserGeometry['Water Level'].x,
                                         UserGeometry['Water Level'].y)
                                   )

#-----------------------------------------------------------------------------
# PLAXIS Negative Interface (for FF boundaries)
#-----------------------------------------------------------------------------
g_i.gotostages()
g_i.activate(g_i.GroundwaterFlowBCs, g_i.InitialPhase)


#-----------------------------------------------------------------------------
# PLAXIS Line Displacements / Motions
#-----------------------------------------------------------------------------
g_i.gotostructures()

ybase = UserGeometry['Bottom'].values
if len(UserGeometry['Bottom'].x) > 2:
    print('Base must be defined by a horizontal line with two points.')        # Can revisit for general procedure
                                                                               # See docs for multiline displ implementation
line_ymin = g_i.line(*zip(ybase[:,0], ybase[:,1]))[-1]
g_i.linedispl(line_ymin)

g_i.LineDisplacements[0].Displacement_x = "Prescribed"   
g_i.LineDisplacements[0].Displacement_y = "Fixed"

if BaseBC == 'fixed':
    g_i.LineDisplacements[0].ux_start = 1.0 # Within motion (expect up + down)
else:
    g_i.LineDisplacements[0].ux_start = 0.5 # Outcropping motion (expect up)

GroundMotion = [None]*N_motion
for count in range(N_motion):

    GroundMotion[count] = g_i.displmultiplier()
    line_ymin.LineDisplacement.LineDisplacement.Multiplierx = GroundMotion[count]
    GroundMotion[count].Signal = "Table"
    GroundMotion[count].Table.set(*Motion[count].itertuples(index=False,name=None))
    GroundMotion[count].DataType = "Accelerations"
    GroundMotion[count].ScalingValue = 9.81
    GroundMotion[count].transform()

PostShake_BaseBC = g_i.displmultiplier()
PostShake_BaseBC.DataType = "Velocities"
PostShake_BaseBC.Frequency = 1

#-----------------------------------------------------------------------------
# PLAXIS FE Mesh
#-----------------------------------------------------------------------------
g_i.gotomesh()

for poly in (g_i.Polygons[:]):
    poly.CoarsenessFactor = 0.15
    
if g_i.SoilPolygons[:]:
    for soilpolygon in (g_i.SoilPolygons[:]):
        soilpolygon.CoarsenessFactor = 0.15 
    
g_i.mesh((0.04),False)
#-----------------------------------------------------------------------------
# PLAXIS Monitoring / Curve Points
#-----------------------------------------------------------------------------
g_i.selectmeshpoints()

for index, row in UserGeometry['Monitoring'].iterrows():
    g_o.addcurvepoint("Node",(row["x"],row["y"]))

g_o.update()

#-----------------------------------------------------------------------------
# PLAXIS Calculation
#-----------------------------------------------------------------------------
# Stage 1 - Ko
g_i.gotostages()  
g_i.InitialPhase.DeformCalcType = "K0 Procedure"
g_i.InitialPhase.PorePresCalcType = "Steady state groundwater flow"
g_i.deactivate(g_i.Dynamics, g_i.InitialPhase)
g_i.InitialPhase.Deform.IgnoreSuction = False #Suction off

# Stage 2 - Initial Stress, Gravity Loading
phase_initial_stress.DeformCalcType = "Plastic"
phase_initial_stress.Identification = "Initial Stress - Gravity Loading"
#phase_initial_stress.PorePresCalcType = "Steady state groundwater flow"       # Can revisit for general procedure
phase_initial_stress.PorePresCalcType = "Phreatic"
g_i.Deformations.BoundaryYMin[phase_initial_stress] = "Vertically fixed"
g_i.deactivate(g_i.Dynamics, phase_initial_stress)
phase_initial_stress.Deform.IgnoreSuction = False #Suction off

# Dynamic and Postliq Stages
#[phase.Identification.value for phase in g_i.Phases]

# Set dynamic stage calculation settings
for i, (dyn_phase) in enumerate(g_i.Phases[dyn_branch_start:dyn_branch_start + N_motion]):
    
    g_i.DynLineDisplacement_1_1.Multiplierx[dyn_phase] = GroundMotion[i] #Displ Mult.
    g_i.DynLineDisplacement_1_2.Multiplierx[dyn_phase] = GroundMotion[i]
    g_i.DynLineDisplacement_1_3.Multiplierx[dyn_phase] = GroundMotion[i]
    
    
    dyn_phase.DeformCalcType = "Dynamic with Consolidation"
    dyn_phase.Identification = MotionNames[i]
    dyn_phase.Deform.UseCavitationCutOff = True # Cavitation cut-off
    dyn_phase.Deform.IgnoreSuction = False      # Ignore suction off
    dyn_phase.ShouldCalculate = True
    dyn_phase.MaxCores = 1
    dyn_phase.Deform.UseDefaultIterationParams = False
    dyn_phase.Deform.TimeStepDetermType = "Manual"

    if len(Motion[i].index) <= 10000:
        dyn_phase.Deform.MaxSteps = len(Motion[i].index)
    else:
        dyn_phase.Deform.SubSteps = 5
        dyn_phase.Deform.MaxSteps = int(np.ceil(len(Motion[i].index)/5))

    dyn_phase.Deform.ResetDisplacementsToZero = True
    dyn_phase.Deform.TimeIntervalSeconds = Motion[i]["Time (s)"].max()
    g_i.Deformations.BoundaryXMin[dyn_phase] = "free"
    g_i.Deformations.BoundaryXMax[dyn_phase] = "free"
    g_i.Deformations.BoundaryYMin[dyn_phase] = "Vertically fixed"
    g_i.activate(g_i.Dynamics,dyn_phase)
    g_i.Dynamics.BoundaryXMin[dyn_phase] = "free-field"
    g_i.Dynamics.BoundaryXMax[dyn_phase] = "free-field"
    if  BaseBC == 'fixed':
        g_i.Dynamics.BoundaryYMin[dyn_phase] = "none"      # Within motion
    else:
        g_i.Dynamics.BoundaryYMin[dyn_phase] = "compliantbase" # Outcropping motion
       
    g_i.activate(g_i.LineDisplacements,dyn_phase)
    g_i.activate(g_i.DynLineDisplacement_1[:],dyn_phase)




# Set post-liq reconsol. stage calculation settings
for i, (post_phase) in enumerate(g_i.Phases[post_branch_start:post_branch_start + N_motion]):
    
    g_i.DynLineDisplacement_1_1.Multiplierx[post_phase] = PostShake_BaseBC
    g_i.DynLineDisplacement_1_2.Multiplierx[post_phase] = PostShake_BaseBC
    g_i.DynLineDisplacement_1_3.Multiplierx[post_phase] = PostShake_BaseBC
    
    
    post_phase.DeformCalcType = "Dynamic with Consolidation"
    post_phase.Identification = MotionNames[i]
    post_phase.Deform.UseCavitationCutOff = True # Cavitation cut-off
    post_phase.Deform.IgnoreSuction = False      # Ignore suction off
    post_phase.ShouldCalculate = True
    post_phase.MaxCores = 1
    post_phase.Deform.UseDefaultIterationParams = False
    post_phase.Deform.TimeStepDetermType = "Manual"
    dyn_phase.Deform.SubSteps = 1
    dyn_phase.Deform.MaxSteps = 1000
    post_phase.Deform.ResetDisplacementsToZero = False
    post_phase.Deform.TimeIntervalSeconds = 10
    g_i.Deformations.BoundaryXMin[post_phase] = "free"
    g_i.Deformations.BoundaryXMax[post_phase] = "free"
    g_i.Deformations.BoundaryYMin[post_phase] = "Vertically fixed"
    g_i.activate(g_i.Dynamics,post_phase)
    g_i.Dynamics.BoundaryXMin[post_phase] = "free-field"
    g_i.Dynamics.BoundaryXMax[post_phase] = "free-field"
    g_i.Dynamics.BoundaryYMin[dyn_phase] = "none"
    g_i.activate(g_i.LineDisplacements,dyn_phase)
    g_i.activate(g_i.DynLineDisplacement_1[:],dyn_phase)



# ...
g_i.gotostages()
g_i.save(r'Model Saves\BHP_EW_Section_FF_PreDyn.p2dx')
#g_i.calculate() 
#g_i.save(r'Model Saves\BHP_EW_Section_FF_Dynamic.p2dx')






# Remaining Items to address:
# Elastic Elements*************** 
#       Maybe expand random field simulation, add option for returning
#       elastic equivalent of intended PM4Sand/Silt...
#       maybe just return the closest at the sides?
# Alpha - Beta
