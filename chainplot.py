# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:22:07 2024

@author: Andoni Agirre Arabolaza
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import warnings
import itertools


def generate_coordinates(center_pos, angle, shape_type, connection_type, a=1):
    """
    Generates two-dimensional shapes in the specified positions and angles, as well as generating the snake
    order of the MPS. Returns the coordinates to generate the polygons and the coordinates of the MPS.

    Parameters
    ----------
    center_pos : Tuple
        Coordinates specifying the coordinate of the center of the shape.
    angle : Float
        Specifies the orientation of the polygon. 
        For hexagons, an angle of 0 has a vertex on top
        For pentagons, an angle of 0 has a vertex on top
    shape_type : string
        Specifies the type of shape: "h" for hexagon, "p" for pentagon.
    connection_type : string
        The order in which the MPS chain transverses the shape.
        Current supported connections are "h", "hup", "hdn", "pup", "pdn"
    a : Float, optional
        The length of the sides of the shape. The default is 1.

    Raises
    ------
    NameError
        Invalid connection_type was specified.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if shape_type == "hexagon" or shape_type == "hex" or shape_type == "h":
        shape_type = "h"
        
        coords = [
            (-np.sqrt(3)*a/2, 0.5),
            (-np.sqrt(3)*a/2, -0.5),
            (0,-a),
            (np.sqrt(3)*a/2,-0.5),
            (np.sqrt(3)*a/2,0.5),
            (0,a)
            ]
        
    elif shape_type == "pentagon" or shape_type == "pent" or shape_type == "p":
        shape_type = "p"
        
        coords = [
            (-0.25*a*(1+np.sqrt(5)), a*np.sqrt((5-np.sqrt(5))/(40))),
            (-0.5*a, -a*np.sqrt((5+2*np.sqrt(5))/(20))),
            (0.5*a, -a*np.sqrt((5+2*np.sqrt(5))/(20))),
            (0.25*a*(1+np.sqrt(5)), a*np.sqrt((5-np.sqrt(5))/(40))),
            (0, a*np.sqrt((5+np.sqrt(5))/(10)))
            ]


    #rotate around the origin first, and only then translate it to the right spot
    def rotate(point, angle):
        angle = angle/180 * np.pi
        x = np.cos(angle)*point[0]-np.sin(angle)*point[1]
        y = np.sin(angle)*point[0]+np.cos(angle)*point[1]
        return (x,y)
    def translate(point, vector):
        return (point[0]+vector[0], point[1]+vector[1])
    def permute(list,permutation):
        if len(list)!=len(permutation):
            warnings.warn("Possible shape missmatch; wrong snake order likely introduced")
        return [list[i] for i in permutation]

    
    coords_rotated = [rotate(vertex, angle=angle) for vertex in coords]
    coords_final = [translate(vertex, vector=center_pos) for vertex in coords_rotated]

    #We also require the coordinates ordered in the mps order!
    #For example the unrotated hexagon centered at the origin, has (in the "h" orientation) snake order
    #coords = [(-np.sqrt(3)*a/2, 0.5),(-np.sqrt(3)*a/2, -0.5), (0,a), (0,-a),(np.sqrt(3)*a/2,0.5),(np.sqrt(3)*a/2,-0.5)]
    h_permutation = [0,1,5,2,4,3]
    hup_permutation = [0,1,2,3,5,4]
    hdn_permutation = [0,1,5,4,3,2]
    
    pdn_permutation = [0,1,4,3,2] #NOTE THIS CHANGES WHICH ONE THE PENTAGON TIP IS, BETTER TO ADJUST SO THAT IT IS ALWAYS (0,a)!
    pup_permutation = [0,1,2,4,3]
    
    #reorder the coordinates to MPS snake order
    if connection_type =="h":
            ordering = h_permutation
    elif connection_type =="hup":
            ordering = hup_permutation
    elif connection_type =="hdn":
        ordering = hdn_permutation
    elif connection_type == "pdn":
        ordering = pdn_permutation
    elif connection_type == "pup":
        ordering = pup_permutation
    elif connection_type == "straight":
        ordering = h_permutation
    else:
        raise NameError('No valid snake order introduced!')
        #print("No valid snake order introduced")
    

    coords_mps = permute(coords_final, ordering)

    
    return coords_final, coords_mps
    

def get_next_center_pos(previous_center, previous_angle, attach_where, type1, type2, a=1):
    """
    Gives the coordinates of the center of the next shape on the chain from the previous
    with the specified shape types and connection types

    Parameters
    ----------
    previous_center : tuple of Floats
        Coordinates of the previous center.
    previous_angle : Float
        Orientation of the previous shape.
    attach_where : String
        Specifies which edge the new shape whould be attached to.
    type1 : String
        Specifies the type of the previous shape.
    type2 : String
        Specifies the type of the next shape.
    a : Float, optional
        The length of the sides of the shapes. The default is 1.

    Raises
    ------
    NameError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    dist = a/2*(np.sqrt(3)+1/(np.tan(np.deg2rad(36))))
    
    if type1=="h" and type2=="p":
        extra_rotation = np.deg2rad(60)
        if attach_where=="below":
            extra_rotation*=-1
        
        
    elif type1=="p" and type2=="h":
        if attach_where == "above":
            extra_rotation = np.deg2rad(54)
            
        elif attach_where == "below":
            extra_rotation = -np.deg2rad(18)
        else:
            raise NameError('Invalid attach location')
    else:
        raise NameError('The specified shapes cannot be attached with the given type combination')
    
    extra_rotation += np.deg2rad(previous_angle)
    x = previous_center[0]+dist*np.cos(extra_rotation)
    y = previous_center[1]+dist*np.sin(extra_rotation)
    return(x,y)


def generate_hp_chain(chaintype=1, reps=1, closinghex=True, extra_rot = 10, a=1):
    """
    Generates the chain. The returned lists have no generated polygons and the MPS
    overcounts the sites shared between shapes.
    Use shapes_to_polys() and clean_mps() to account for those.

    Parameters
    ----------
    chaintype : Int, optional
        Only chain type 1 is supported now. The default is 1.
    reps : Int, optional
        Number of times the shapes are repeated. The default is 1.
    closinghex : Bool, optional
        Specifies whether a closing hexagon should be added to the end of the chain. The default is True.
    extra_rot : Float, optional
        Rotates the entire chain by the specified amount in deg (Not implemented yet). The default is 10.
    a : Float, optional
        Length of the edges of the shapes. The default is 1.

    Returns
    -------
    shapes : list
        Coordinates of all vertices shapes of the shapes.
    mps_list : list
        Coordinates of .

    """
    if chaintype == 1:
        angles= list(np.tile([0, 42, 42-18, 42-18-78], reps))
        centers= [(0,0)]
        polys = list(np.tile(["h", "p", "h", "p"], reps))

        attach = list(np.tile(["above", "below", "below", "above"], reps))
        connections = list(np.tile(["hup", "pdn", "hdn", "pup"], reps))
        
        if closinghex:
            polys.append("h")
            angles.append(0)
            connections.append("hup")

    elif chaintype == 2:
        angles = list(np.tile([0, -18, -36, -54], reps))
        centers= [(0,0)]
        polys = list(np.tile(["h", "p", "h", "p"], reps))
        attach = list(np.tile(["straight", "below", "straight", "below"], reps))
        connections = list(np.tile(["h", "pdn", "h", "pdn"], reps))
        
        if closinghex:
            polys.append("h")
            angles.append(0)
            connections.append("h")
        
        
        
    else:
        print("Unsupported chain type")
    
    shapes = []
    mps_list = []
    
    shape, mps = generate_coordinates(centers[0], angles[0], polys[0], connection_type=connections[0], a=a);
    shapes.append(shape)
    mps_list.append(mps)
    
    for i in range(len(angles)-1):
        centers.append(get_next_center_pos(centers[i], angles[i], attach[i], polys[i], polys[i+1]))
        shape, mps = generate_coordinates(centers[i+1], angles[i+1], polys[i+1], connection_type=connections[i+1], a=a);
        shapes.append(shape)
        mps_list.append(mps)
        
    return shapes, mps_list


def clean_mps(mps_list, flatten=True, to_plot=True):
    """
    Removes the duplicate coordinates of the MPS which represnt the same site

    Parameters
    ----------
    mps_list : list
        List of mps coordinates.
    flatten : Bool, optional
        Specifies whether the final list should be flattened or the shape-structure preserved. The default is True.

    Returns
    -------
    new_mps : list
        List of mps coordinates with no duplicates.

    """
    new_mps = []
    num_shapes = len(mps_list)
    last_reached=False
    for i, shape in enumerate(mps_list):
        if i != (num_shapes - 1):
            assert(last_reached == False), "MPS error, shape iterated over after last shape"
            new_mps.append(shape[:-2]) #remove last two (repeated) sites
        else:
            new_mps.append(shape)
            last_reached = True
            
    if flatten:
        new_mps = [x for xs in new_mps for x in xs]
    
    if to_plot:
        x=[i[0] for i in new_mps]
        y=[i[1] for i in new_mps]
        new_mps = [x,y]
    return new_mps


def shapes_to_polys(shapes, kwargs={'fc':'white', 'ec':'k'}):
    #Generates the polygons from the shape lists
    #Use sx.add_patch() to add these to axis ax
    polys = []
    for s in shapes:
        polys.append(Polygon(s, **kwargs))
    return polys

    

def generate_chain(chaintype=1, reps=1, closinghex=True, extra_rot = 10, a=1, flatten=True, kwargs={'fc':'white', 'ec':'k'}):
    #Generates the full ladder chain of type chaintype from scratch
    #Use ax.add_patch(poly) to add each element of polys to axis ax
    #Use ax.plot(*mps_clean) to plot the mps
    #Any observables can be plotted straightforwardly if they follow the mps order
    sh, mp = generate_hp_chain(chaintype=chaintype, reps=reps, closinghex=closinghex, extra_rot = extra_rot, a=a)
    polys = shapes_to_polys(sh, kwargs=kwargs)
    mps_clean = clean_mps(mp, flatten=flatten, to_plot=True)
    return polys, mps_clean


def get_midpoint(a, b):
    x = 0.5*(a[0]+b[0])
    y = 0.5*(a[1]+b[1])
    return (x,y)

def generate_edge_coordinates(mps, a):
    midpoints = []
    for combination in itertools.combinations(mps, 2):
        p1 = combination[0]
        p2 = combination[1]
        if np.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2)<(1.01*a):
            midpoints.append(get_midpoint(p1, p2))
    return midpoints

def generate_mps_bond_coordinates(mps, site1, site2):
    return get_midpoint(mps[site1], mps[site2])

def generate_lattice_bonds(geometry):
    #Takes something like ["h", "pup"]
    
    bonds_hex = [(1,2), 2, 2, 2]
    bonds_pentdn = [(1,2), 3, 1] #down means the next hex goes below, ie, the free vertex is up!
    bonds_pentup = [(1,3), 1, 2]
    bonds_hexdn = [(1,2), 4, 1, 1]
    bonds_hexup = [(1,4), 1, 1, 2]
    bonds_close = [1] #These always close with 1, unless system is periodic, in that case we should work with a size n-1 unit cell! And end with (-(N-2), -(N-2))

    
    def symtobond(sym):
        if sym == "h":
            return bonds_hex
        elif sym == "pup":
            return bonds_pentup
        elif sym == "pdn":
            return bonds_pentdn
        elif sym == "c":
            return bonds_close  
        elif sym == "hup":
            return bonds_hexup  
        elif sym == "hdn":
            return bonds_hexdn  


    """
    periodic_bonds = []
    

    N_shapes = length(geometry)
    N_sites_period =  
    for char in geometry
        if char == "h"
            append!(periodic_bonds, bonds_hex)
        else if char == "p↑"
            append!()
        elseif char

    end
    """

    #bonds = [symtobond(sym); for sym in geometry]
    #GENERATE (NON-PERIODIC) BONDS
    bonds = []
    for g in geometry:
        print("g:", g)
        bonds = [*bonds, *symtobond(g)]
        print(bonds)
    

    if geometry[-1] != "c":
        bonds = [*bonds, *bonds_close]
    
    #@assert length(bonds)==N_sites Not true, it will have at least one less!
    return bonds

def generate_lattice_from_bonds(bonds):
    lat = []
    bond_counter = 1
    
    def makeiterable(x):
        if type(x) == int:
            return list([x])
        else:
            return x
    for i in range(len(bonds)):
        for j in makeiterable(bonds[i]):
            print("Adding bond from site i=", i, " to i+", j,"=",i+j)
            #j tells us how many sites ahead we should connect to from i
            #for example, for a single hexagon in ring-snake configuration we would have bonds = [(1,2),(2),(2),(2),(1)]
            lat.append((i, i+j)) #no coordinates for now. sites_list[i]..., sites_list[i+2]..., "")
            #println(lat[bond_counter])
            bond_counter+=1
            print(bond_counter)

    return lat



"""
geometry = ["hup", "pdn", "hdn", "pup", "hup"]
a = generate_lattice_bonds(geometry)
b = generate_lattice_from_bonds(a)
#Use this by calling chain() and then ax.add_patch(polys[:]).
"""
