#"bread and butter" calculations
import numpy as np
from Data_Processing import *
from Data_Managment import *

[Airfoil_Database, Airfoil_Database_strings] = get_airfoil_database()

body_list_path_base = 'saves/'
all_mass_path_base = 'saves/'
balance_list_base_path = 'saves/'

#finds all unknown FS, BL, and WL based on section span and sweep angle
def calculate_planform_unknowns(Planform_Geo, Body_name):

    num_sections = len(Planform_Geo)

    new_array = [0,0,0,0,0,0,0,0]

    for s in range(num_sections):

        Input_Type = Planform_Geo[s][7]

        if s < num_sections-1:
            next_Input_Type = Planform_Geo[s+1][7]

        #If user wants to manually input all data (shouo
        if Input_Type == 0:

            if s < num_sections-1:
                if next_Input_Type == 0:

                    FS_LE = Planform_Geo[s][0]
                    BL = Planform_Geo[s][1]
                    WL = Planform_Geo[s][2]
                    span = Planform_Geo[s][3]

                    next_FS_LE = Planform_Geo[s+1][0]
                    next_WL = Planform_Geo[s+1][2]
                    LE_Sweep = rad2deg(np.arctan((next_FS_LE-FS_LE)/span))
                    Dihedral = rad2deg(np.arctan((next_WL-WL)/span))

                    Hinge_index = Planform_Geo[s][6]

                else:

                    FS_LE = Planform_Geo[s][0]
                    BL = Planform_Geo[s][1]
                    WL = Planform_Geo[s][2]
                    span = Planform_Geo[s][3]

                    LE_Sweep = Planform_Geo[s][4]
                    Dihedral = Planform_Geo[s][5]

                    Hinge_index = Planform_Geo[s][6]
            else:

                FS_LE = Planform_Geo[s][0]
                BL = Planform_Geo[s][1]
                WL = Planform_Geo[s][2]
                span = Planform_Geo[s][3]

                LE_Sweep = Planform_Geo[s][4]
                Dihedral = Planform_Geo[s][5]

                Hinge_index = Planform_Geo[s][6]

        #calculate FS, BL, WL from previous sweep and dihedral values
        elif Input_Type == 1:

            # Pull in planform data and calculate it if values are not explicit

            last_FS_LE = Planform_Geo[s-1][0]
            FS_LE = last_FS_LE + np.tan(deg2rad(Planform_Geo[s - 1][4])) * Planform_Geo[s - 1][3]
            Planform_Geo[s][0] = FS_LE

            last_BL = Planform_Geo[s-1][1]
            BL = last_BL + Planform_Geo[s - 1][3]
            Planform_Geo[s][1] = BL

            last_WL = Planform_Geo[s-1][2]
            WL = last_WL + np.tan(deg2rad(Planform_Geo[s - 1][5])) * Planform_Geo[s - 1][3]
            Planform_Geo[s][2] = WL

            span = Planform_Geo[s][3]

            if next_Input_Type == 1:

                LE_Sweep = Planform_Geo[s][4]
                Dihedral = Planform_Geo[s][5]

            else:

                next_FS_LE = Planform_Geo[s + 1][0]
                next_WL = Planform_Geo[s + 1][2]
                LE_Sweep = rad2deg(np.arctan((next_FS_LE - FS_LE) / span))
                Dihedral = rad2deg(np.arctan((next_WL - WL) / span))

            Hinge_index = Planform_Geo[s][6]

        new_row = [FS_LE,BL,WL,span,LE_Sweep,Dihedral,Hinge_index,Input_Type]

        new_array = np.vstack((new_array, new_row))

    new_array = new_array[1:] #chop off the top row

    return new_array

def Calculate_Surface(aero_Body):

    planform_Geo = aero_Body.planform_Geo
    section_Geo = aero_Body.section_Geo

    Body_name = aero_Body.name
    body_mirror = aero_Body.mirror

    num_sections = len(section_Geo)
    num_surfaces = len(planform_Geo) - 1

    Total_Area = 0
    Total_Span = 0
    MAC_area = 0
    MACy_area = 0
    MACx_area = 0

    Area_mat = [0, 0, 0]
    Hinge_mat = [0, 0, 0, 0]
    for s in range(num_sections - 1):

        # Pull in section data
        chord = section_Geo[s][0]
        chord_next = section_Geo[s + 1][0]
        twist = section_Geo[s][1]
        hinge_x = section_Geo[s][2]
        hinge_x_next = section_Geo[s + 1][2]
        hinge_z = section_Geo[s][3]
        hinge_r = section_Geo[s][4]

        big_chord = chord + chord_next

        #eventually pull ac data from airfoil database
        ac = 0.25
        ac_next = 0.25

        FS_LE = planform_Geo[s][0]
        FS_LE_next = planform_Geo[s + 1][0]
        BL = planform_Geo[s][1]
        WL = planform_Geo[s][2]
        span = planform_Geo[s][3]

        ac_FS = ((ac * chord) + FS_LE)
        ac_FS_next = ((ac_next * chord_next) + FS_LE_next)

        ac_FS_xchord = ac_FS * chord
        ac_FS_next_xchord = ac_FS_next * chord_next

        Total_Span += span

        hinge_index = planform_Geo[s][6]  # which hinge index this hinge is attributed to (-1 means it is not a true hinge)

        Surface_Area = ((chord + chord_next) / 2) * span

        Total_Area += Surface_Area #add to the total surface area of the body

        Centroid_span = (span * ((2 * chord_next) + chord)) / (3 * (chord_next + chord))
        Centroid_BL = BL + Centroid_span

        #find the FS of the aerodynamic center based on trapazoid using the centroid span found above
        ac_FS_x = ac_FS+ ((ac_FS_next - ac_FS) / span) * (Centroid_span)

        ac_FS_area = ac_FS_x * Surface_Area
        MACx_area += ac_FS_area

        Centroid_BL_area = Centroid_BL * Surface_Area
        MACy_area += Centroid_BL_area

        mac = chord + ((Centroid_span / span) * (chord_next - chord)) #find mean aerodynamic chord of section
        mac_area = mac * Surface_Area #create a factor that we can scale in order to find the final MAC

        MAC_area += mac_area

        hinge_chordwise1 = (1 - hinge_x) * chord #chordwise size of the hinge
        hinge_chordwise2 = (1 - hinge_x_next) * chord_next #chordwise size of next section

        hinge_spanwise = span

        hinge_area = ((hinge_chordwise1 + hinge_chordwise2) / 2) * hinge_spanwise

        hinge_centroid_span = (span * ((2 * hinge_chordwise2) + hinge_chordwise1)) / (3 * (hinge_chordwise2 + hinge_chordwise1))
        hinge_centroid_BL = hinge_centroid_span + BL

        hinge1_C_FS = FS_LE + (chord*(hinge_x + ((1-hinge_x)/2)))
        hinge2_C_FS = FS_LE_next + (chord_next*(hinge_x_next + ((1-hinge_x_next)/2)))

        hinge_centroid_FS = (hinge1_C_FS+hinge2_C_FS)/2


        hinge_list = [hinge_index,hinge_area,hinge_centroid_FS,hinge_centroid_BL]

        data_list = [Surface_Area, hinge_area, hinge_index]

        Area_mat = np.vstack((Area_mat, data_list))
        Hinge_mat = np.vstack((Hinge_mat, hinge_list))

    Area_mat = Area_mat[1:] #[section area, hinge area, hinge index]
    Hinge_mat = Hinge_mat[1:]

    Hinge_Array = Boil_Hinges(Hinge_mat, aero_Body)

    MAC = MAC_area / Total_Area
    MACy = MACy_area / Total_Area
    MACx = MACx_area / Total_Area
    MAC_ac = 0.25

    body_mirror = np.double(body_mirror)
    if body_mirror == 1:
        Total_Area = Total_Area * 2

    MAC_vec = [MACx,MACy,MAC,MAC_ac,Total_Area]

    return [MAC_vec, Area_mat, Hinge_Array]

def Boil_Hinges(hinge_mat, aero_Body):

    num_hinges = aero_Body.num_hinge

    unique_hinges = np.unique(hinge_mat[:, 0])
    unique_areas = np.zeros(num_hinges+1)
    unique_FS = np.zeros(num_hinges+1)
    unique_BL = np.zeros(num_hinges + 1)
    for n in range(len(hinge_mat)):
        hinge_sec_area = hinge_mat[n][1]
        hinge_index = int(hinge_mat[n][0])
        hinge_C_FS = hinge_mat[n][2]
        hinge_C_BL = hinge_mat[n][3]

        C_FS_A = hinge_C_FS * hinge_sec_area
        C_BL_A = hinge_C_BL * hinge_sec_area

        unique_areas[hinge_index + 1] = unique_areas[hinge_index + 1] + hinge_sec_area
        unique_FS[hinge_index + 1] = unique_areas[hinge_index + 1] + C_FS_A
        unique_BL[hinge_index + 1] = unique_areas[hinge_index + 1] + C_BL_A

    for i in range(num_hinges+1):
        hinge_area = unique_areas[i]
        hinge_FS = unique_FS[i]
        hinge_BL = unique_BL[i]

        if hinge_area != 0:
            unique_FS[i] = hinge_FS / hinge_area
            unique_BL[i] = hinge_BL / hinge_area
        else:
            unique_FS[i] = -999
            unique_BL[i] = -999

    hinge_array = np.vstack((unique_areas,unique_FS,unique_BL)).transpose()

    return hinge_array

#Calculate the CG of all of the points, also calculate the
def calculate_CG(save_name):

    all_mass_path = '{0}{1}/Mass_Inertial_Properties.csv'.format(all_mass_path_base, save_name)

    all_mass = process_string_csv(all_mass_path,0,0)
    mass_intertial_properties = np.double(all_mass[:,1:])

    num_items = len(all_mass)
    total_mass = sum(mass_intertial_properties[:,0])

    sum_CG_FS = 0
    sum_CG_BL = 0
    sum_CG_WL = 0

    sum_Ixx = 0
    sum_Iyy = 0
    sum_Izz = 0
    sum_Iyz = 0
    sum_Ixz = 0
    sum_Ixy = 0

    for n in range(num_items):

        row = mass_intertial_properties[n]

        mass = row[0]
        item_FS = row[1]
        item_BL = row[2]
        item_WL = row[3]

        sum_CG_FS += mass * item_FS
        sum_CG_BL += mass * item_BL
        sum_CG_WL += mass * item_WL

    CG_FS = sum_CG_FS / total_mass
    CG_BL = sum_CG_BL / total_mass
    CG_WL = sum_CG_WL / total_mass

    for i in range(num_items):

        row = mass_intertial_properties[i]

        mass = row[0]
        item_FS = row[1]
        item_BL = row[2]
        item_WL = row[3]

        sum_Ixx += mass * ((item_BL - CG_BL)**2 + (item_WL - CG_WL)**2)
        sum_Iyy += mass * ((item_FS - CG_FS)**2 + (item_WL - CG_WL)**2)
        sum_Izz += mass * ((item_FS - CG_FS)**2 + (item_BL - CG_BL)**2)
        sum_Iyz += mass * ((item_BL - CG_BL) * (item_WL - CG_WL))
        sum_Ixz += mass * ((item_FS - CG_FS) * (item_WL - CG_WL))
        sum_Ixy += mass * ((item_FS - CG_FS) * (item_BL - CG_BL))

    inertial_vec = [total_mass, CG_FS, CG_BL, CG_WL, sum_Ixx, sum_Iyy, sum_Izz, sum_Iyz, sum_Ixz, sum_Ixy]

    return inertial_vec

#calulate the AC of the aircraft based on the areas and ac's of each horizontal body
def calculate_aircraft_ac(save_name):

    balance_list_path = '{0}{0}/Aero_Balance_List.csv'.format(balance_list_base_path, save_name)
    body_list_path = '{0}{1}/Body_List.csv'.format(body_list_path_base, save_name)

    balance_list = process_csv(balance_list_path, 0, 0)
    body_list = process_string_csv(body_list_path, 0, 0)
    orientation_list = np.double(body_list[:,1])
    mirror_list = np.double(body_list[:,2])

    area_list = balance_list[:,4]

    sum_area = 0
    sum_ac_FS_area = 0
    sum_ac_BL_area = 0
    for n in range (len(orientation_list)):

        orientation = orientation_list[n]
        mirror = mirror_list[n]
        area = area_list[n]
        ac_FS = balance_list[n][0]
        ac_BL = balance_list[n][1]

        if orientation == 0:

            sum_area += area
            sum_ac_FS_area += ac_FS * area

            if mirror == 1:
                sum_ac_BL_area = sum_ac_BL_area
            else:
                sum_ac_BL_area += ac_BL * area

    sum_area = round(sum_area, 2)
    total_ac_FS = round(sum_ac_FS_area / sum_area, 2)
    total_ac_BL = round(sum_ac_BL_area / sum_area, 2)

    return [sum_area, total_ac_FS, total_ac_BL]
