import pickle
import numpy as np
import pandas as pd
from shapely.geometry import (MultiLineString, LineString, MultiPolygon, 
                              Polygon, GeometryCollection, Point)
from shapely.ops import split
import matplotlib.pylab as plt
from scipy import interpolate


def read_excel_poly(filepath):
    ## Packs excel file polyzone inputs conviently into a dictionary
    
    df = pd.read_excel(filepath,
                       engine='openpyxl',
                       skiprows=0,
                       )
    Zones = []
    Names = []
    TableIDs = [] 

    table_loc = 1
    header_loc = 2
    Ncol = 2 
    N = int(len(df.columns))
    for i in range(0, N, Ncol):
        
        df_zone = df.iloc[:,i:i+Ncol].copy().dropna()
        name = df.iloc[:,i:i+Ncol].columns.values[1]
        table_id = df.iloc[:,i:i+Ncol].loc[0].values[1]
        
        if ('Unnamed' not in name):
            new_columns = df_zone.iloc[table_loc,:].values
            df_zone = df_zone.iloc[header_loc:,:]
            df_zone.columns = new_columns
            
            Zones.append(df_zone)
            Names.append(name)
            TableIDs.append(table_id)
            
    poly_dict = dict(zip(Names,Zones))
            
    return poly_dict, (Zones, Names, TableIDs)



def mesh_grid(UserGeometry, GenerateZones, NxNy, case='grid', plot_mesh=True):
    """
    ## Generates collection of mesh grid polyzones for PLAXIS instantiation
    ## from user-specified geometry
    
    Parameters
    ----------
    UserGeometry : DICT
        Dictionary of polyzone coordinates (pandas.core.frame.DataFrame). 
        Should use keys compatible with GenerateZones
    GenerateZones : LIST
        List of units to generate. Should use keys compatible with UserGeometry
    NxNy : TUPLE
        Number of columns and rows (Nx, Ny).
    case : STR, optional
       'grid' for regular square gridding, or 'topo' to follow topographic
       boundaries. The default is 'grid'.
    plot_mesh : TYPE, optional
        Plot zoning. The default is True.
    Returns
    -------
    PolyMeshCollections : Dict
        Dictionary of polyzone collections
        (shapely.geometry.multipolygon.MultiPolygon) with keys compatible with
        UserGeometry and Generate Zones.
    PolyCoordCollections : Dict
        Dictionary of polyzone coordinates
        (shapely.geometry.multipolygon.MultiPolygon) with keys compatible with
        UserGeometry and Generate Zones.
    """

    
    if case == 'grid':  # Create regular square gridding
        

        PolyMeshCollections = {}
        PolyCoordCollections = {}
        for zone in GenerateZones:
            
            points = list(zip(UserGeometry[zone].x, 
                              UserGeometry[zone].y)
                          )
            
            nx, ny = NxNy[zone]     # number of columns and rows
            
            polygon = Polygon(points)
            
            minx, miny, maxx, maxy = polygon.bounds
            dx = (maxx - minx) / nx  # width of a small part
            dy = (maxy - miny) / ny  # height of a small part
            horizontal_splitters = [LineString([(minx, miny + i*dy),
                                                (maxx, miny + i*dy)]) for i in range(ny)]
            vertical_splitters = [LineString([(minx + i*dx, miny),
                                              (minx + i*dx, maxy)]) for i in range(nx)]
            splitters = horizontal_splitters + vertical_splitters
            
            
            results = polygon
            for splitter in splitters:
                results = MultiPolygon(split(results, splitter))
            
               
            #parts = [list(part.exterior.coords) for part in result.geoms]
            #print(parts)
            coords = [list(poly.exterior.coords) for poly in results]
            
            PolyMeshCollections[zone] = results
            PolyCoordCollections[zone] = coords
            
            if plot_mesh:
                for poly in results:
                    plt.plot(*poly.exterior.xy)
                    
    elif case == 'topo':  # Create gridding which follows topographic boundaries
        
        PolyMeshCollections = {}
        PolyCoordCollections = {}
        for zone in GenerateZones:
            
            points = list(zip(UserGeometry[zone].x, 
                              UserGeometry[zone].y)
                          )
            
            nx, ny = NxNy[zone]     # number of columns and rows
            
            polygon = Polygon(points)
            
            minx, miny, maxx, maxy = polygon.bounds
            dx = (maxx - minx) / nx  # width of a small part
            
            # Get surface topo vertices
            x_topo_vertices = np.unique(np.sort(polygon.exterior.xy[0]))
            
            # Create multilines between top/bottom boundaries
            horizontal_lines = []
            Xs = []
            Ys = []
            for j in range(0,ny+1):
                xs = []
                ys = []
                for i in range(0,len(x_topo_vertices)):
                    
                    # Create vertical linestring, and find intersections w/ top and
                    # bottom of polygon for getting element height
                    v_line = LineString([(x_topo_vertices[i], miny),
                                          (x_topo_vertices[i], maxy)]
                                        )
                    #print( polygon.intersection(v_line) )
                    
                    xy_intersec = polygon.intersection(v_line).xy
                    ybott = min(xy_intersec[1])
                    ytop = max(xy_intersec[1])
                    dy = abs(ytop - ybott)/ny         # height of a small part
                
                    x = x_topo_vertices[i]
                    y = ytop - j*dy
                    #print(f"x:{x}, y:{y}, ybott:{ybott}, ytop:{ytop}, dy:{dy}")
                    
                    xs.append(x)
                    ys.append(y)
                
                if len(xs) > 1:
                    LS = LineString(zip(xs,ys))
                    horizontal_lines.append(LS)
                    
                Xs += xs
                Ys += ys
            
            # Create polygons from linestrings
            Polys = []
            for i in range(len(horizontal_lines)-1):
                l1 = horizontal_lines[i].coords
                l2 = horizontal_lines[i+1].coords
                Polys.append( Polygon([*list(l1), *list(l2)[::-1]]))
            #parts = [list(poly.exterior.coords) for poly in Polys]
                
            
            # Vertical splitters
            vertical_splitters = [LineString([(minx + i*dx, miny),
                                              (minx + i*dx, maxy)]) for i in range(nx)]
            splitters = vertical_splitters
            
            # Apply vertical splitters to each horizontal layer polygon
            results = []
            for poly in Polys:
                result = poly
                for splitter in splitters:
                    result = MultiPolygon(split(result, splitter))
                results.append(result)
                
            # flatten the list
            flattened_polys = []
            for i in range(len(results)):
                for j in range(len(results[i])):
                    flattened_polys.append(results[i][j])
            results = MultiPolygon(flattened_polys)
               
            #parts = [list(part.exterior.coords) for part in results.geoms]
            #print(parts)
            coords = [list(poly.exterior.coords) for poly in results]
            
            PolyMeshCollections[zone] = results #GeometryCollection(results)
            PolyCoordCollections[zone] = coords
        
        if plot_mesh:
            for result in results:
                for poly in result:
                    plt.plot(*poly.exterior.xy)
                    plt.plot(*polygon.exterior.xy)
            plt.plot(Xs,Ys,'o', ms=2.5)
        
    return PolyMeshCollections, PolyCoordCollections



def map_DR(load_RF_path, realiz_idx, start_elev,
            PolyMeshCollections, PolyCoordCollections, zone, 
            DR_min = 0.35, DR_max = 0.85,
            preview=False):
    ## Function for loading pickled DR random field and interpolating
    ## values (LEGACY--to be updated with random field generation files)
    
    # Load pickle file
    with open(load_RF_path,'rb') as f:
        out_dict = pickle.load(f)
    
    # Define the interpolating function
    # .. grid and DR
    grid = out_dict['GridStack'].copy()
    randomfield = out_dict['RFFinal'][realiz_idx].copy()
    # .. replace nan values with zero for interpolation
    randomfield[np.isnan(randomfield)] = 0 
    x = np.arange(grid.x1['start'],
                  grid.x1['start'] + 
                  grid.x1['cell dim']*grid.x1['cell size'],
                  grid.x1['cell size'])
    y = np.arange(grid.x2['start'],
                  grid.x2['start'] + 
                  grid.x2['cell dim']*grid.x2['cell size'],
                  grid.x2['cell size'])
    # .. correct to elevation
    y = start_elev - y
    z = randomfield
    # .. create interpolating fcn
    f_DR = interpolate.interp2d(x, y, z, kind='cubic')
    
    # .. interpolate within all model polygon zones using points inside polygon    # Refine for sampling more pts 
    N = len(PolyCoordCollections[zone])
    xc, yc, z_new = np.zeros(N), np.zeros(N), np.zeros(N)
    for i, (poly, xys) in enumerate(zip(PolyMeshCollections[zone],
                                        PolyCoordCollections[zone])
                                    ):
        
        xc[i] = poly.centroid.x
        yc[i] = poly.centroid.y
        
        xys += [(xc[i], yc[i])]
        
        zi = np.array([f_DR(xy[0],xy[1]) for xy in xys])
        zi = zi[zi>0].mean() # Average positive values
        z_new[i] = zi
    
    # Apply relative density bounds
    z_new[z_new < DR_min] = DR_min 
    z_new[z_new > DR_max] = DR_max     

    # Preview
    if preview:
        fig = plt.figure(figsize=(20,5))
        plt.plot(xc,yc,'ko',ms=0.5) # gridpoints
        plt.tricontourf(xc,yc,z_new) # use for irregular points, 3 1d arrays
        plt.show()

    return z_new, (xc, yc), f_DR


def get_zone_xyc(multipolygoncollection):
    ## Return representative x,y coordinates of elements for stress calcs, 
    ## either the centroid or shapely representative point
    
    N = len(multipolygoncollection)
    xc, yc = np.zeros(N), np.zeros(N)
    for i, poly in enumerate(multipolygoncollection):
        
        smaller = Point((poly.centroid.x, poly.centroid.y))
        if poly.contains(smaller):
            # Ideal for dense meshing/rectangular elements, when elements
            # are widely spaced or irregular, centroid can fall outside 
            xc[i] = poly.centroid.x
            yc[i] = poly.centroid.y   
        
        else:
            # Garanteed to fall inside polygon
            xc[i] = poly.representative_point().x
            yc[i] = poly.representative_point().y 
            
    return xc, yc


def get_RGB_number(R, G, B):
    ## Convert RGB values to PLAXIS compatible colors
    
    # get colour number from RGB using BIT LEFT SHIFT
    iB = B<<16 # left shift 16 bits for Blue
    iG = G<<8  # left shift  8 bits for Green
    iR = R     # left shift  0 bits for Red
    return iB + iG + iR