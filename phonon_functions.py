def get_delta_angle(structure, range_angle, conv_matrix=[]):
    # Return a list of delta angles from a range of %

    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    lattice = np.matmul(conv_matrix,structure.lattice.matrix)

    angle = Structure(lattice,[],[]).lattice.alpha

    delta_angle = []
    for i in np.arange(range_angle[0],range_angle[1]+range_angle[2]/2,range_angle[2]):
        delta_angle.append(angle*i)

    return delta_angle

def get_delta_lattice(structure, range_lattice, conv_matrix=[]):
    # Return a list of delta lattice from a range of %

    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    lattice = np.matmul(conv_matrix,structure.lattice.matrix)

    lattice = Structure(lattice,[],[]).lattice.a

    delta_lattice = []
    for i in np.arange(range_lattice[0],range_lattice[1]+range_lattice[2]/2,range_lattice[2]):
        delta_lattice.append(lattice*i)
        print(lattice,i)

    return delta_lattice

def set_angle(structure, delta_angle, conv_matrix=[]):
    # structure: primitive cell structure
    # delta_angle: delta on the conventional cell angle
    # conv_matrix: conversion matrix from the primitive to conventional cell
    
    import numpy as np
    from pymatgen.core.structure import Structure

    to_rad = np.pi/180
    to_deg = 180/np.pi


    if len(conv_matrix) = 0:
        conv_matrix = np.identity(3)

    conv_matrix_inv = np.linalg.inv(conv_matrix)
    lattice = np.matmul(conv_matrix,structure.lattice.matrix)


    A = lattice[0,:]
    B = lattice[1,:]
    C = lattice[2,:]

    len_A = np.linalg.norm(lattice[0,:])
    len_B = np.linalg.norm(lattice[1,:])
    len_C = np.linalg.norm(lattice[2,:])

    A_norm = A/len_A
    B_norm = B/len_B
    C_norm = C/len_C

    alpha = np.arccos(np.dot(B_norm,C_norm))*(180/np.pi)
    beta = np.arccos(np.dot(A_norm,C_norm))*(180/np.pi)
    gamma = np.arccos(np.dot(A_norm,B_norm))*(180/np.pi)

    vector_sum = A_norm + B_norm + C_norm
    len_vector_sum = np.linalg.norm(vector_sum)

    vector_sum_norm = vector_sum/len_vector_sum

    D_a = np.cross(np.cross(A_norm,vector_sum_norm),A_norm)
    D_b = np.cross(np.cross(B_norm,vector_sum_norm),B_norm)
    D_c = np.cross(np.cross(C_norm,vector_sum_norm),C_norm)

    len_D_a = np.linalg.norm(D_a)
    len_D_b = np.linalg.norm(D_b)
    len_D_c = np.linalg.norm(D_c)

    D_a_norm = D_a/len_D_a
    D_b_norm = D_b/len_D_b
    D_c_norm = D_c/len_D_c

    alpha_t = np.arccos(np.dot(A_norm,vector_sum_norm))*(180/np.pi)
    beta_t = np.arccos(np.dot(B_norm,vector_sum_norm))*(180/np.pi)
    gamma_t = np.arccos(np.dot(C_norm,vector_sum_norm))*(180/np.pi)

    len_vector_sum_p = np.sqrt(3 + 2*(np.cos((alpha+delta_angle)*to_rad) +
                                        np.cos((beta+delta_angle)*to_rad) +
                                        np.cos((gamma+delta_angle)*to_rad)))


    alpha_p = (np.arccos((  1+ np.cos((beta+delta_angle)*to_rad) + np.cos((gamma+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg
    beta_p = (np.arccos((  1+ np.cos((alpha+delta_angle)*to_rad) + np.cos((gamma+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg
    gamma_p = (np.arccos((  1+ np.cos((alpha+delta_angle)*to_rad) + np.cos((beta+delta_angle)*to_rad) ) /len_vector_sum_p))*to_deg


    delta_alpha_t = (alpha_t-alpha_p)*to_rad
    delta_beta_t = (beta_t-beta_p)*to_rad
    delta_gamma_t = (gamma_t-gamma_p)*to_rad

    #print(alpha,alpha_t,alpha_p,delta_alpha_t)

    A1 = (np.cos(delta_alpha_t)*A_norm + np.sin(delta_alpha_t)*D_a_norm) #* len_A
    B1 = (np.cos(delta_beta_t)*B_norm + np.sin(delta_beta_t)*D_b_norm) #* len_B
    C1 = (np.cos(delta_gamma_t)*C_norm + np.sin(delta_gamma_t)*D_c_norm) #* len_C

    new_lattice_conv = np.array([A1,B1,C1])

    new_lattice_conv[0,:] = new_lattice_conv[0,:] *(len_A)
    new_lattice_conv[1,:] = new_lattice_conv[1,:] *(len_B)
    new_lattice_conv[2,:] = new_lattice_conv[2,:] *(len_C)

    new_lattice_prim = np.matmul(conv_matrix_inv,new_lattice_conv)


    return Structure(new_lattice_prim,structure.atomic_numbers,structure.frac_coords)  

def set_cell(structure, delta_cell, conv_matrix=[]):
    # structure: primitive cell structure
    # delta_cell: delta on the conventional cell vectors
    # conv_matrix: conversion matrix from the primitive to conventional cell
    
    import numpy as np
    from pymatgen.core.structure import Structure

    if len(conv_matrix) == 0:
        conv_matrix = np.identity(3)

    conv_matrix_inv = np.linalg.inv(conv_matrix)
    lattice = np.matmul(conv_matrix,structure.lattice.matrix)


    A = lattice[0,:]
    B = lattice[1,:]
    C = lattice[2,:]

    len_A = np.linalg.norm(lattice[0,:])
    len_B = np.linalg.norm(lattice[1,:])
    len_C = np.linalg.norm(lattice[2,:])


    ###
    lattice_new = np.zeros((3,3))
    lattice_new[0,:] = (lattice[0,:]/len_A) * (len_A+delta_cell)
    lattice_new[1,:] = (lattice[1,:]/len_B) * (len_B+delta_cell)
    lattice_new[2,:] = (lattice[2,:]/len_C) * (len_C+delta_cell)
    
    prim_lattice_new = np.matmul(conv_matrix_inv,lattice_new)
    
    return Structure(prim_lattice_new,structure.atomic_numbers,structure.frac_coords)

  