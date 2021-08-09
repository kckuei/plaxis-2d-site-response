# import os
#import pickle
#import numpy as np
#import pandas as pd
# from shapely.geometry import (MultiLineString, LineString, MultiPolygon, 
#                               Polygon, GeometryCollection)
# from shapely.ops import split
import matplotlib.pylab as plt
#from scipy import interpolate

import utils as ut



#-----------------------------------------------------------------------------
# Unit Tests
#-----------------------------------------------------------------------------


def unit_test_grid():
    ## Example - generate polyzone mesh with regular grid
    
    UserGeometry, _ = ut.read_excel_poly(r"User Inputs\User_Geometry.xlsx") # Load polyzones
    
    
    GenerateZones = ['Sand', 'Clay', 'Till']        # Specify zones to generate
    
    NxNy = {'Sand':(40, 20),                # Specify number of columns and rows
            'Clay':(5, 7), 
            'Till':(1, 1),
            }    
    
    PolyMeshCollections, PolyCoordCollections = ut.mesh_grid(UserGeometry,
                                                            GenerateZones,
                                                            NxNy,
                                                            case='grid',
                                                            plot_mesh=False)
    
    for zone in PolyMeshCollections.keys():
        results = PolyMeshCollections[zone]
        for poly in results:
            plt.plot(*poly.exterior.xy)
    plt.show()
            
    for zone in PolyCoordCollections.keys():
        result = PolyCoordCollections[zone]
        for coords in result:
            list_of_tuples = coords
            x, y = map(list, zip(*list_of_tuples))
            plt.plot(x,y)
    plt.show()
    return None
    
def unit_test_topo():
    ## Example - generate polyzone mesh following topographic bounds
    
    UserGeometry, _ = ut.read_excel_poly(r"User Inputs\User_Geometry.xlsx") # Load polyzones
    
    
    GenerateZones = ['Sand', 'Clay', 'Till']        # Specify zones to generate
    
    NxNy = {'Sand':(40, 20),                # Specify number of columns and rows
            'Clay':(5, 7), 
            'Till':(1, 1),
            }

    PolyMeshCollections, PolyCoordCollections = ut.mesh_grid(UserGeometry,
                                                            GenerateZones,
                                                            NxNy,
                                                            case='topo',
                                                            plot_mesh=False)
    
    for zone in PolyMeshCollections.keys():
        results = PolyMeshCollections[zone]
        for poly in results:
            plt.plot(*poly.exterior.xy)
    #plt.axis('equal')
    plt.show()
            
    for zone in PolyCoordCollections.keys():
        result = PolyCoordCollections[zone]
        for coords in result:
            list_of_tuples = coords
            x, y = map(list, zip(*list_of_tuples))
            plt.plot(x,y)
    plt.show()
    return None

def centroid_comparison_A():
    ## Shows spatial comparison of choice of using 'representative points', 
    ## average of polygon bounds, and centroid for initial stress calcs
    ## Centroid ideal in this example
    
    UserGeometry, _ = ut.read_excel_poly(r"User Inputs\User_Geometry.xlsx") # Load polyzones
    
    
    GenerateZones = ['Sand', 'Clay', 'Till']        # Specify zones to generate
    
    NxNy = {'Sand':(40, 20),                # Specify number of columns and rows
            'Clay':(5, 7), 
            'Till':(1, 1),
            }

    PolyMeshCollections, PolyCoordCollections = ut.mesh_grid(UserGeometry,
                                                            GenerateZones,
                                                            NxNy,
                                                            case='topo',
                                                            plot_mesh=False)
    
    # Plot comparison
    zone = 'Sand'
    fig, ax = plt.subplots(figsize = (12,8))
    for poly in PolyMeshCollections[zone]:
            
            minx, miny, maxx, maxy = poly.bounds
            
            plt.plot((minx+maxx)/2, (miny+maxy)/2,
                     's', ms = 7.5,
                     )
            plt.plot(poly.representative_point().x,
                     poly.representative_point().y,
                     '^', ms = 7.5,
                     )
            plt.plot(poly.centroid.x, poly.centroid.y,
                     'o', ms = 7.5,
                     )
            plt.plot(*poly.exterior.xy,lw=0.5)
    
    plt.plot([],[],'ks', mfc='w', ms = 7.5,label = 'Average Bounds'
             )
    plt.plot([],[],'k^', mfc='w', ms = 7.5,label = 'Representative'
             )
    plt.plot([],[],'ko', mfc='w', ms = 7.5,label = 'Centroid'
             )
    plt.xlim(1500,2000)
    plt.ylim(-70,10 )
    plt.legend(loc='upper right')
    return None

def centroid_comparison_B():
    ## Shows spatial comparison of choice of using 'representative points', 
    ## average of polygon bounds, and centroid for initial stress calcs
    ## Representative ideal in this example
    
    UserGeometry, _ = ut.read_excel_poly(r"User Inputs\User_Geometry.xlsx") # Load polyzones
    
    
    GenerateZones = ['Sand', 'Clay', 'Till']        # Specify zones to generate
    
    NxNy = {'Sand':(5, 20),                # Specify number of columns and rows
            'Clay':(5, 7), 
            'Till':(1, 1),
            }

    PolyMeshCollections, PolyCoordCollections = ut.mesh_grid(UserGeometry,
                                                            GenerateZones,
                                                            NxNy,
                                                            case='topo',
                                                            plot_mesh=False)
    
    # Plot comparison
    zone = 'Sand'
    fig, ax = plt.subplots(figsize = (12,8))
    for poly in PolyMeshCollections[zone]:
            
            minx, miny, maxx, maxy = poly.bounds
            
            plt.plot((minx+maxx)/2, (miny+maxy)/2,
                     's', ms = 7.5,
                     )
            plt.plot(poly.representative_point().x,
                     poly.representative_point().y,
                     '^', ms = 7.5,
                     )
            plt.plot(poly.centroid.x, poly.centroid.y,
                     'o', ms = 7.5,
                     )
            plt.plot(*poly.exterior.xy,lw=0.5)
    
    plt.plot([],[],'ks', mfc='w', ms = 7.5,label = 'Average Bounds'
             )
    plt.plot([],[],'k^', mfc='w', ms = 7.5,label = 'Representative'
             )
    plt.plot([],[],'ko', mfc='w', ms = 7.5,label = 'Centroid'
             )
    plt.xlim(1500,2000)
    plt.ylim(-70,10 )
    plt.legend(loc='upper right')
    return None
    
unit_test_grid()
unit_test_topo()
centroid_comparison_A()
centroid_comparison_B()