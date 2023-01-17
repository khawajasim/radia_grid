#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:22:22 2022

@author: khawajasim
"""
import numpy
import numpy as np
import json 
import mercantile
import matplotlib.pyplot as plt
from shapely.geometry import mapping # box, shape, mapping, Point, Polygon

def plot_grid_polygon(grid_polygon):
    fig, ax = plt.subplots()
    fig.set_size_inches(15, 15)
    
    
    for arc in grid_polygon:
        coords = arc.exterior.xy
        ax.plot(coords[0], coords[1]) #color = 'black'
        
    fig.tight_layout()
    
    return ax


def plot_stress(stress_data, which_data = 'CFS'):
    """
    Plots stress as heatmap.
    Stress_data: 2D numpy.array()
                [dx,dy, depth, stress1, stress2, earthquakes, ....]
    """
    fig, ax = plt.subplots()

    # im = ax.imshow(data, **kwargs)
    dx = numpy.sort(numpy.unique(stress_data[:,0]))
    data_2d = []
    for xx in dx:
        data_xx = stress_data[stress_data[:,0] == xx]
        data_xx = data_xx[data_xx[:,1].argsort()]
        if which_data == 'CFS':
            data_2d.append(data_xx[:,3])
        if which_data == 'static':
            data_2d.append(data_xx[:,4])
    data_2d = numpy.array(data_2d)
    #plt.imshow(data_2d, cmap='Blues')
    im = ax.imshow( numpy.log10(data_2d), cmap='Reds' )
    #plot negative data. Because log can be plotted in negative.
    data_2d[data_2d>0] = 0
    data_2d = abs(data_2d)
    im = ax.imshow(numpy.log10(data_2d), cmap='Blues')
        # tick_locs = dx[ax.get_yticks().astype(int)[1:]]
        # ax.set_xticklabels([-100, -75, -50, -25, 0, 25, 50, 75, 100])
    ax.set_yticks([])
    ax.set_xticks([])
      
    return im

# (stress_data, Num_eqs=[137, 0], for_all_positive=False, for_all_negative=False)
def generate_active_cells(stress_data, Num_eqs=[137, 0], for_all_positive=False, for_all_negative=False):
    """
    Generate earthquakes in a cell in such a way that exactly "Num_eqs" number of cells should be active.
    stress_data: n x 5
                [x, y, depth, coulomb_stress, reference_stress]
    Num_eqs = [num earthquakes at positive, num earthquakes at negative]
    """
    #separate reference stress based on positive and negative coulomb stress
    stress_positive = stress_data[stress_data[:,3]>0]
    stress_negative = stress_data[stress_data[:,3]<0]

    if for_all_positive: 
        #This is an extreme case, where we assign earthquake to all positive stress
        events_positive_stress = numpy.ones(len(stress_positive))
        events_negative_stress = numpy.zeros(len(stress_negative))
    
    elif for_all_negative:
        #This is an extreme case, where we assign earthquake to all negative stress
        events_positive_stress = numpy.zeros(len(stress_positive))
        events_negative_stress = numpy.ones(len(stress_negative))
            
    else:
        #Generate GIVEN NUMBER of positive and negative stress (whichever required)
        events_positive_stress = numpy.zeros(len(stress_positive))
        events_negative_stress = numpy.zeros(len(stress_negative))
    
        #----Generate earthquakes for positive stress using Reference model (Distance based stress)
        eqs=0
        probmap_positive = stress_positive[:,4]
        probmap_positive = numpy.cumsum(probmap_positive) / sum(probmap_positive)
        while eqs < Num_eqs[0]:
            random_num = numpy.random.uniform(0,1)
            loc = numpy.searchsorted(probmap_positive, random_num)
            if events_positive_stress[loc] == 0:
                numpy.add.at(events_positive_stress, loc, 1)
                eqs = eqs+1
        print('---Positive stress earthquakes generated: ', eqs)
#        
   
        #----Generate earthquakes for Negative stress using Reference model (Distance based stress)
        eqs = 0
        probmap_negative = stress_negative[:,4]
        probmap_negative = numpy.cumsum(probmap_negative) / sum(probmap_negative)
        while eqs < Num_eqs[1]:
            random_num = numpy.random.uniform(0,1)
            loc = numpy.searchsorted(probmap_negative, random_num)
            if events_negative_stress[loc] == 0:
                numpy.add.at(events_negative_stress, loc, 1)
                eqs = eqs+1
        print('---Negative stress earthquakes generated: ', eqs)
#        stress_negative_eqs = numpy.column_stack((stress_negative, events_negative_stress)) 

    #Appending all positive and negative earthquakes to respective stress data.
    stress_positive_eqs = numpy.column_stack((stress_positive, events_positive_stress))
    stress_negative_eqs = numpy.column_stack((stress_negative, events_negative_stress))
    stress_eqs = numpy.row_stack((stress_positive_eqs, stress_negative_eqs))
    
    return stress_eqs


def write_geojson_feature(grid_polygons,feature, filename):
    """
    Write coordsniates along with feature
    This functions takes in Shape.Geometery type list of shapes (polygon or point) and writes
    them in the format of geojason file. The geojason file will be then used to plot in QGIS

    Parameters
    ----------
    shape_list : List of Shapely.geometery
        Bounding boxes or points or lines that are used for display in QGIS
    filename : string
        Intended name of GEOJSON file

    Returns 
    -------
    GEOJSON file

    """
    shape_list=[] #Converting to the proper format. 
    for poly in grid_polygons:
        shape_list.append(mapping(poly))
        
        
    geo_shape = {"type":"FeatureCollection",
                  "features":[]}
    for bb, ff in zip(shape_list, feature):
        feature = {"type":"Feature",
              "properties":{"feature1": ff,
                            },
                  "geometry":bb};
        geo_shape['features'].append(feature)

#   poly['features'][0]['geometry']['coordinates']

    with open(filename+'.geojson', 'a') as f:
        json.dump(geo_shape, f)
    
    return geo_shape



def best_mcc_f1_threshold(mcc_f1_data):
    """
    nx3: ndarrray
        [mcc, f1, threshold]
        
    """
    dist = numpy.sqrt((mcc_f1_data[:,0]-1)**2 + (mcc_f1_data[:,1]-1)**2)
    return mcc_f1_data[numpy.argmin(dist)][2]


def trapz_interval(y, x=None, dx=1.0, axis=-1):
    """
    Calculate area under curve using Trapezoidal rule. 
    Calculate mean area, area considering left bounds, area consdiering right bouds
    inputs: 
        y axis data (TPR in our case)
        x axis data (FPR in our case)
    """
    if x is None:
        d = dx
    else:
#        x = asanyarray(x)
        if x.ndim == 1:
            d = numpy.diff(x)
            # reshape to correct shape
            shape = [1]*y.ndim
            shape[axis] = d.shape[0]
            d = d.reshape(shape)
        else:
            d = numpy.diff(x, axis=axis)
    nd = y.ndim
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1, None)
    slice2[axis] = slice(None, -1)
    area = (d * (y[tuple(slice1)] + y[tuple(slice2)]) / 2.0).sum(axis)
    left_bound_area = (d * y[tuple(slice2)]).sum(axis)
    right_bound_area = (d * y[tuple(slice1)]).sum(axis)
    return area, left_bound_area, right_bound_area

def calc_F1_R_SenSepc(stress, earthquakes):
    F1_value= []
    SenSpec_value = []
    Rscore_value = []
    for CFSth in numpy.sort(stress):
#        print(CFSth)
        TP = len(stress[((stress>=CFSth) & (earthquakes > 0))])
        FP = len(stress[((stress>=CFSth) & (earthquakes == 0))])
        FN = len(stress[((stress <CFSth) & (earthquakes > 0))])
        TN = len(stress[((stress <CFSth) & (earthquakes == 0))])
#        print('TP :',TP)
#        print('FP :' ,FP)
#        print('FN :',FN)
#        print('TN :',TN)
        if FN > 0 and TN > 0:
#            MCC = (((TP*TN)-(FP*FN))/((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))**0.5 +1)/2 #MCC rescaled between 0-1
#            print(MCC)
            F1_value.append((2*TP)/((2*TP)+FP+FN))       #F1-measure
            SenSpec_value.append((TN/(TN+FP)) * (TP/(TP+FN)))  #specificity*Sensitivty (Not working)
            Rscore_value.append((((TP*TN) - (FP*FN)) / ((TP+FN)*(FP+TN)) +1)/2) #R-score
#            MCC_F1.append([F1, MCC, CFSth])    #[x-axis data, y-axis data, threshold]    
    
    return F1_value, SenSpec_value, Rscore_value


def _sample_earthquakes(N, prob, grid_bounds):
    
    """
    Use probablity map "prob" to simulate "N" earthquakes location (x,y) for each cell characterized by grid_bounds:[lon1,lat1,lon2,lat2]
    as  bottom_left and top_right corners.
    
    grid_bounds: nx[lon1,lat1,lon2,lat2]
    """
    x = np.sort(np.random.uniform(0, 1, N))
    prob = prob/sum(prob)
    pcum = np.cumsum(prob)
    Ncum = len(pcum)
#    print("Ncum : ",Ncum)
    xloc = np.zeros(N)
    yloc = np.zeros(N)
    j = 0
    for i in range(N):
        while j+1<Ncum and pcum[j+1]<=x[i]:
            j += 1

        xloc[i] = np.random.uniform(grid_bounds[j,0], grid_bounds[j,2])
        yloc[i] = np.random.uniform(grid_bounds[j,1], grid_bounds[j,3])  #---
    return numpy.column_stack((xloc, yloc))


def generate_earthquakes_xy(grid_bounds,  stress_data, Num_eqs=[137, 0]):
    """
    Simulate earthquakes locations (x,y) for each cell, separately for Positive and Negative Coulom stress.
    The coordinates (x,y) of simulated earthquakes are generated for each cell characterized by xbounds and ybounds as origin
    coordinates with cells size of dh. 
    Inputs:
        grid_bounds: nx4
                [lon1,lat1,lon2,lat2]
        stress_data: n x 2
                [x, y, depth, coulomb_stress, reference_stress]
        dh: cell size to compute bounds.
        Num_eqs = [num earthquakes at positive, num earthquakes at negative]
    Outputs:
        nx[x,y] 
        
    """
    #separate reference stress based on positive and negative coulomb stress
    stress_positive = stress_data[stress_data[:,0]>0]
    grid_bounds_positive = grid_bounds[stress_data[:,0]>0]
    
    stress_negative = stress_data[stress_data[:,0]<0]
    grid_bounds_negative = stress_data[stress_data[:,0]<0]
    #Generate GIVEN NUMBER of positive and negative stress (whichever required)
    # events_positive_stress = numpy.zeros(len(stress_positive))
    # events_negative_stress = numpy.zeros(len(stress_negative))
    
    #----Generate earthquakes for positive stress using Reference model (Distance based stress)
    
    probmap_positive = stress_positive[:,1]
    probmap_positive = probmap_positive / sum(probmap_positive)
    # _sample_earthquakes(N, prob, grid_bounds)
    stress_positive_eqs = _sample_earthquakes(Num_eqs[0], probmap_positive, grid_bounds_positive)
    
    #----Generate earthquakes for Negative stress using Reference model (Distance based stress)
    probmap_negative = stress_negative[:,2]
    probmap_negative = probmap_negative / sum(probmap_negative)
    stress_negative_eqs = _sample_earthquakes(Num_eqs[1], probmap_negative, grid_bounds_negative)
     
    #Appending all positive and negative earthquakes to respective stress data.
   
    return numpy.row_stack((stress_positive_eqs, stress_negative_eqs))



def generate_gridded_eqs(grid, stress_data, Num_eqs=[137, 0], model_type= 'ref', for_all_positive=False, for_all_negative=False):
    """
    Generate earthquakes in a cell in such a way that exactly "Num_eqs" number of cells should be active.
    Inputs:
        grid: shapely.polygons corresponding to the grid cells.
        stress_data: n x 2
                [coulomb_stress, reference_stress]
        Num_eqs = [num earthquakes at positive, num earthquakes at negative]
        model_type = 'ref' or 'coulomb': to generate the synthetics.
    Outputs:
        grid: numpy.array(shapely.polygons) 
                Arranged as positive and negative.
        stress_eqs: nx3
                [coulomb_stress, reference_stress, eqs]
                arranged as positive, negative stress.
        
    """
    #separate reference stress based on positive and negative coulomb stress
    grid = numpy.array(grid)
    stress_positive = stress_data[stress_data[:,0]>0]
    grid_positive = grid[stress_data[:,0]>0]
    stress_negative = stress_data[stress_data[:,0]<0]
    grid_negative = grid[stress_data[:,0]<0]
    if for_all_positive: 
        #This is an extreme case, where we assign earthquake to all positive stress
        events_positive_stress = numpy.ones(len(stress_positive))
        events_negative_stress = numpy.zeros(len(stress_negative))
    
    elif for_all_negative:
        #This is an extreme case, where we assign earthquake to all negative stress
        events_positive_stress = numpy.zeros(len(stress_positive))
        events_negative_stress = numpy.ones(len(stress_negative))
            
    else:
        #Generate GIVEN NUMBER of positive and negative stress (whichever required)
        events_positive_stress = numpy.zeros(len(stress_positive))
        events_negative_stress = numpy.zeros(len(stress_negative))
        
        
        # ****Choose Stress type to generate synthetics.
        if model_type == 'ref':
            probmap_positive = stress_positive[:,1]
            probmap_negative = stress_negative[:,1]
        elif model_type == 'coulomb':
            probmap_positive = stress_positive[:,0]
            probmap_negative = stress_negative[:,0]
        else: 
            print('Choose the correct model type for catalog generation')
            quit()
            
        
        #----Generate earthquakes for positive stress using Reference model (Distance based stress)
        probmap_positive = numpy.cumsum(probmap_positive) / sum(probmap_positive)    
        
        random_numbers = numpy.random.rand(Num_eqs[0])
        # find insertion points using binary search inserting to satisfy a[i-1] <= v < a[i]
        pnts = numpy.searchsorted(probmap_positive, random_numbers, side='right')
        # create simulated catalog by adding to the original locations
        numpy.add.at(events_positive_stress, pnts, 1)   
        
        
        
        #----Generate earthquakes for Negative stress using Reference model (Distance based stress)
        probmap_negative = numpy.cumsum(probmap_negative) / sum(probmap_negative)
        
        random_numbers = numpy.random.rand(Num_eqs[1])
        # find insertion points using binary search inserting to satisfy a[i-1] <= v < a[i]
        pnts = numpy.searchsorted(probmap_negative, random_numbers, side='right')
        # create simulated catalog by adding to the original locations
        numpy.add.at(events_negative_stress, pnts, 1)   
        
#        stress_negative_eqs = numpy.column_stack((stress_negative, events_negative_stress)) 

    #Appending all positive and negative earthquakes to respective stress data.
    stress_positive_eqs = numpy.column_stack((stress_positive, events_positive_stress))
    stress_negative_eqs = numpy.column_stack((stress_negative, events_negative_stress))
    stress_eqs = numpy.row_stack((stress_positive_eqs, stress_negative_eqs))
    grid = numpy.append(grid_positive, grid_negative)
    
    return grid, stress_eqs

def grid_origin_coordinates(qk):
    """
    Parameters
    ----------
    qk : Array of Strings
        Quad keys.

    Returns
    -------
    grid_coords : Array of floats
        Origin Coordinates of Grid formed by Quad keys boxes
        [Longitude Latitude]

    """
    
    grid_lat =[]
    grid_lon = []
    for i in range(len(qk)):
        grid_lon = numpy.append(grid_lon,mercantile.bounds(mercantile.quadkey_to_tile(qk[i])).west)
        grid_lat = numpy.append(grid_lat,mercantile.bounds(mercantile.quadkey_to_tile(qk[i])).south)

    grid_coords = numpy.column_stack([grid_lon, grid_lat])
    return grid_coords

def grid_top_right_coordinates(qk):
    """
    Parameters
    ----------
    qk : Array of Strings
        Quad keys.

    Returns
    -------
    grid_top_coords : Array of floats
        Top-Right Coordinates of Grid formed by Quad keys boxes
        [Longitude Latitude]

    """
    
    grid_lat =[]
    grid_lon = []
    for i in range(len(qk)):
        grid_lon = numpy.append(grid_lon,mercantile.bounds(mercantile.quadkey_to_tile(qk[i])).east)
        grid_lat = numpy.append(grid_lat,mercantile.bounds(mercantile.quadkey_to_tile(qk[i])).north)

    grid_top_coords = numpy.column_stack([grid_lon, grid_lat])
    return grid_top_coords
