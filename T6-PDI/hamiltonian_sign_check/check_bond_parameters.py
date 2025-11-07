import numpy as np
from rotate_onto_xy import *

file = 'T6-single-crystal.xyz'
file_xy='T6-single-crystal-xy-plane.xyz'

def read_array(file):

    with open(file) as geom_file:
        geom_list = geom_file.readlines()
        geom_list = geom_list[2:]

        mass_list = []
        for element in geom_list:
            
            if element[0] == 'S': mass_list.append(32)
            elif element[0] == 'H': mass_list.append(1)
            elif element[0] == 'O': mass_list.append(16)
            elif element[0] == 'N': mass_list.append(14)
            else: mass_list.append(12)

        new_geom_list = []
        for element in geom_list:
            new_geom_list.append(element[1:].split())

        geom_array = np.array(new_geom_list)

        vector_float = np.vectorize(float)
        geom_array = vector_float(geom_array) 

    return np.array(geom_array), np.array(mass_list)

def calculate_angles(coordinates):
    num_atoms = len(coordinates)
    angles = np.zeros((num_atoms, num_atoms, num_atoms))
    for i in range(num_atoms):
        for j in range(num_atoms):
            for k in range(j + 1, num_atoms):
                if i != j and i != k:
                    vector_ij = coordinates[j] - coordinates[i]
                    vector_ik = coordinates[k] - coordinates[i]
                    cos_angle = np.dot(vector_ij, vector_ik) / (np.linalg.norm(vector_ij) * np.linalg.norm(vector_ik))
                    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))  # Clip for numerical stability
                    angles[i, j, k] = np.degrees(angle)
                    angles[i, k, j] = angles[i, j, k]  # The angle is the same regardless of the order
    return angles

def calculate_distances(coordinates):

    number_atoms = len(coordinates)
    distances = np.zeros((number_atoms, number_atoms))

    for index in range(len(coordinates)):
        distance_array = np.zeros(number_atoms)
        
        for index2 in range(len(coordinates)):
            distance = np.linalg.norm(coordinates[index] - coordinates[index2])
            distance_array[index2] = distance

        distances[:, index] = distance_array
    
    return distances



geom, mass_array = read_array(file)
xy_geom = transform_molecule(geom, mass_array)

xy_geom_file, mass_array_file = read_array(file_xy)
print(xy_geom_file - xy_geom)

angles = calculate_angles(geom)
#xy_angles = calculate_angles(xy_geom)

distances = calculate_distances(geom)
#xy_distances = calculate_distances(xy_geom)

#for number in range(len(distances)):
 #   print(distances[:,number] - xy_distances[:,number])
