#All plotting functions
# all return matplotlib 'fig' data types

import matplotlib.pyplot as plt
import numpy as np

from Data_Managment import *
from Base_Functions import *
from Base_Calculations import *

airfoil_database_path = gen_filename('SHARED', 'airfoil_database', 0)

Airfoil_Database = process_csv(airfoil_database_path, 3, 0)
Airfoil_Database_strings = process_string_csv(airfoil_database_path, 0, 0)
Airfoil_Database_strings = Airfoil_Database_strings[0, :]

[Airfoil_Database, Airfoil_Database_strings] = get_airfoil_database()

#options_mat = process_string_csv('config/Options.csv', 0, 0)

#was 'process_plot_body'
# takes data from body and solves for all data needed to make graphs
def process_plot_body(Section_Geo, Planform_Geo, Airfoil_Geo, orientation, output_name):

    #length_unit = options_mat[0][0]
    length_unit = 'ft' #locked for testing

    length_unit_list_nums = [1, 12,1/25.4,1/2.54,1/0.0254]
    length_unit_list = ['in', 'ft', 'mm', 'cm', 'm']

    length_to_in_index = length_unit_list.index(length_unit)
    length_to_in = length_unit_list_nums[length_to_in_index]

    ChordLine_Geo = np.array([[0,0],[1,0]])

    num_sections = len(Section_Geo)
    num_surfaces = len(Planform_Geo) - 1

    all_elements = [['', '', '']]
    LE_Guide = ['', '', '']
    TE_Guide = ['', '', '']
    Plotting_Data_Mat = np.array([0, 0, 0, 0], dtype='object')
    for s in range(num_sections):

        # Pull in section data
        chord = Section_Geo[s][0] * length_to_in
        twist = Section_Geo[s][1]
        hinge_x = Section_Geo[s][2]
        hinge_z = Section_Geo[s][3]
        hinge_r = Section_Geo[s][4]

        # Pull in planform data and calculate it if values are not explicit
        FS_LE = Planform_Geo[s][0] * length_to_in
        BL = Planform_Geo[s][1] * length_to_in
        WL = Planform_Geo[s][2] * length_to_in

        # find the airfoil index and select the airfoil points corresponding to that index
        index = np.where(Airfoil_Database_strings == Airfoil_Geo[s])[0]
        ind_col = index[0], index[0]+1
        Airfoil_Pts = Airfoil_Database[:, ind_col]
        Airfoil_Pts = process_airfoil_pts(Airfoil_Pts)

        # make curve for hinge circle
        hinge_circle = circle_pts(hinge_x, hinge_z, hinge_r)
        hinge_circle = np.multiply(hinge_circle, chord)

        hinge_circle[0] = hinge_circle[0] + FS_LE
        hinge_circle[1] = hinge_circle[1] + WL

        # scale rotate and translate the airfoil points and the chordline
        Final_Airfoil_Pts = rotate_scale_translate(Airfoil_Pts,chord,twist,hinge_x,hinge_z,FS_LE,WL)
        Final_Chordline_Pts = rotate_scale_translate(ChordLine_Geo,chord,twist,hinge_x,hinge_z,FS_LE,WL)

        # put spline and guide points into 3D arrays
        BL_array = np.array([BL] * len(Final_Airfoil_Pts[0]))

        if orientation == '0':
            D_Pts = Final_Airfoil_Pts[0], BL_array, Final_Airfoil_Pts[1]
        elif orientation == '1':
            D_Pts = Final_Airfoil_Pts[0], Final_Airfoil_Pts[1], BL_array

        BL_array_Guides = np.array([BL] * len(Final_Chordline_Pts[0]))
        Guide_Pts = Final_Chordline_Pts[0], BL_array_Guides, Final_Chordline_Pts[1]

        # finalize guide points and add them to the guide list
        LE_guide_Pts = Guide_Pts[0][0], Guide_Pts[1][0], Guide_Pts[2][0]
        TE_guide_Pts = Guide_Pts[0][1], Guide_Pts[1][1], Guide_Pts[2][1]

        if orientation == '0':
            LE_guide_Pts = np.vstack((LE_guide_Pts[0],LE_guide_Pts[1],LE_guide_Pts[2])).T
            TE_guide_Pts = np.vstack((TE_guide_Pts[0], TE_guide_Pts[1], TE_guide_Pts[2])).T
        elif orientation == '1':
            LE_guide_Pts = np.vstack((LE_guide_Pts[0], LE_guide_Pts[2], LE_guide_Pts[1])).T
            TE_guide_Pts = np.vstack((TE_guide_Pts[0], TE_guide_Pts[2], TE_guide_Pts[1])).T

        LE_Guide = np.vstack((LE_Guide,LE_guide_Pts))
        TE_Guide = np.vstack((TE_Guide, TE_guide_Pts))

        new_Plotting_Data = np.array([BL / length_to_in, Final_Airfoil_Pts / length_to_in, Final_Chordline_Pts / length_to_in, hinge_circle / length_to_in], dtype='object')

        Plotting_Data_Mat = np.vstack((Plotting_Data_Mat,new_Plotting_Data))


    Plotting_Data_Mat = Plotting_Data_Mat[1:]
    Plotting_Data_Mat = Plotting_Data_Mat

    # [BL, Airfoil_Pts, Chordline_Pts, Hinge_Circle_Pts]
    return Plotting_Data_Mat

#updated version of 'process_plot_body' that runs off object type
def process_plot_body_object(Aero_Body):

    Planform_Geo = Aero_Body.planform_Geo
    Section_Geo = Aero_Body.section_Geo
    Airfoil_Geo = Aero_Body.airfoil_Geo

    orientation = Aero_Body.orientation


    #length_unit = options_mat[0][0]
    length_unit = 'ft' #locked for testing

    length_unit_list_nums = [1, 12,1/25.4,1/2.54,1/0.0254]
    length_unit_list = ['in', 'ft', 'mm', 'cm', 'm']

    length_to_in_index = length_unit_list.index(length_unit)
    length_to_in = length_unit_list_nums[length_to_in_index]

    ChordLine_Geo = np.array([[0,0],[1,0]])

    num_sections = len(Section_Geo)
    num_surfaces = len(Planform_Geo) - 1

    all_elements = [['', '', '']]
    LE_Guide = ['', '', '']
    TE_Guide = ['', '', '']
    Plotting_Data_Mat = np.array([0, 0, 0, 0], dtype='object')
    for s in range(num_sections):

        # Pull in section data
        chord = Section_Geo[s][0] * length_to_in
        twist = Section_Geo[s][1]
        hinge_x = Section_Geo[s][2]
        hinge_z = Section_Geo[s][3]
        hinge_r = Section_Geo[s][4]

        # Pull in planform data and calculate it if values are not explicit
        FS_LE = Planform_Geo[s][0] * length_to_in
        BL = Planform_Geo[s][1] * length_to_in
        WL = Planform_Geo[s][2] * length_to_in

        # find the airfoil index and select the airfoil points corresponding to that index
        index = np.where(Airfoil_Database_strings == Airfoil_Geo[s])[0]
        ind_col = index[0], index[0]+1
        Airfoil_Pts = Airfoil_Database[:, ind_col]
        Airfoil_Pts = process_airfoil_pts(Airfoil_Pts)

        # make curve for hinge circle
        hinge_circle = circle_pts(hinge_x, hinge_z, hinge_r)
        hinge_circle = np.multiply(hinge_circle, chord)

        hinge_circle[0] = hinge_circle[0] + FS_LE
        hinge_circle[1] = hinge_circle[1] + WL

        # scale rotate and translate the airfoil points and the chordline
        Final_Airfoil_Pts = rotate_scale_translate(Airfoil_Pts,chord,twist,hinge_x,hinge_z,FS_LE,WL)
        Final_Chordline_Pts = rotate_scale_translate(ChordLine_Geo,chord,twist,hinge_x,hinge_z,FS_LE,WL)

        # put spline and guide points into 3D arrays
        BL_array = np.array([BL] * len(Final_Airfoil_Pts[0]))

        if orientation == '0':
            D_Pts = Final_Airfoil_Pts[0], BL_array, Final_Airfoil_Pts[1]
        elif orientation == '1':
            D_Pts = Final_Airfoil_Pts[0], Final_Airfoil_Pts[1], BL_array

        BL_array_Guides = np.array([BL] * len(Final_Chordline_Pts[0]))
        Guide_Pts = Final_Chordline_Pts[0], BL_array_Guides, Final_Chordline_Pts[1]

        # finalize guide points and add them to the guide list
        LE_guide_Pts = Guide_Pts[0][0], Guide_Pts[1][0], Guide_Pts[2][0]
        TE_guide_Pts = Guide_Pts[0][1], Guide_Pts[1][1], Guide_Pts[2][1]

        if orientation == '0':
            LE_guide_Pts = np.vstack((LE_guide_Pts[0],LE_guide_Pts[1],LE_guide_Pts[2])).T
            TE_guide_Pts = np.vstack((TE_guide_Pts[0], TE_guide_Pts[1], TE_guide_Pts[2])).T
        elif orientation == '1':
            LE_guide_Pts = np.vstack((LE_guide_Pts[0], LE_guide_Pts[2], LE_guide_Pts[1])).T
            TE_guide_Pts = np.vstack((TE_guide_Pts[0], TE_guide_Pts[2], TE_guide_Pts[1])).T

        LE_Guide = np.vstack((LE_Guide,LE_guide_Pts))
        TE_Guide = np.vstack((TE_Guide, TE_guide_Pts))

        new_Plotting_Data = np.array([BL / length_to_in, Final_Airfoil_Pts / length_to_in, Final_Chordline_Pts / length_to_in, hinge_circle / length_to_in], dtype='object')

        Plotting_Data_Mat = np.vstack((Plotting_Data_Mat,new_Plotting_Data))


    Plotting_Data_Mat = Plotting_Data_Mat[1:]
    Plotting_Data_Mat = Plotting_Data_Mat

    # [BL, Airfoil_Pts, Chordline_Pts, Hinge_Circle_Pts]
    return Plotting_Data_Mat

def Plot_Airfoil_Section(s, Plotting_Data):

    BL = Plotting_Data[s][0]
    Final_Airfoil_Pts = Plotting_Data[s][1]
    Final_Chordline_Pts = Plotting_Data[s][2]
    hinge_circle = Plotting_Data[s][3]

    plt.figure(figsize=(20, 20))

    plt.plot(Final_Airfoil_Pts[0], Final_Airfoil_Pts[1], color='b', linewidth=1)
    plt.plot(Final_Chordline_Pts[0], Final_Chordline_Pts[1], color='black', linewidth=.5, linestyle='--')
    plt.plot(hinge_circle[0], hinge_circle[1], color='red', linewidth=.5, linestyle='--')

    # format the graph
    plt.axis('equal')

    fig = plt.gcf()
    plt.xlabel('FS (X)')
    plt.ylabel('WL (Z)')

    return fig

def Plot_All_Airfoil_Section(Plotting_Data):

    # start plotting

    num_hinge = len(Plotting_Data)

    color_mat = gen_color_mat((0.1, 0.1, 1), (1, 0.1, 0.25), int(num_hinge))

    plt.figure(figsize=(20, 20))

    for s in range(len(Plotting_Data)):

        # color finder
        section_color = color_mat[int(s)]

        BL = Plotting_Data[s][0]
        Final_Airfoil_Pts = Plotting_Data[s][1]
        Final_Chordline_Pts = Plotting_Data[s][2]
        hinge_circle = Plotting_Data[s][3]

        plt.plot(Final_Airfoil_Pts[0], Final_Airfoil_Pts[1], color=section_color, linewidth=1)
        plt.plot(Final_Chordline_Pts[0], Final_Chordline_Pts[1], color='black', linewidth=.5, linestyle='--')

    # format the graph
    plt.axis('equal')

    fig = plt.gcf()

    plt.xlabel('FS (X)')
    plt.ylabel('WL (Z)')

    return fig

def Plot_Airfoil(points):

    plt.clf()

    plt.plot(points[:,0], points[:,1])

    plt.axis('equal')

    fig = plt.gcf()

    return fig

def Plot_Planform(aero_Body):

    planform_Geo = aero_Body.planform_Geo
    section_Geo = aero_Body.section_Geo
    mirror = aero_Body.mirror

    num_sections = len(section_Geo)
    num_hinge = max(planform_Geo[:, 6]) + 1

    color_mat = gen_color_mat((0.1, 0.1, 1), (1, 0.1, 0.25), int(num_hinge))

    MAC_vec = Calculate_Surface(aero_Body)
    MAC_vec = MAC_vec[0]

    MAC_LE_FS = (MAC_vec[2] * (-MAC_vec[3])) + MAC_vec[0]
    MAC_TE_FS = (MAC_vec[2] * (1 - MAC_vec[3])) + MAC_vec[0]

    mirror = float(mirror)

    plt.figure(figsize=(20, 20)) #force the figure to be large (figure will still scale to window size in GUI)

    #draw aerodynamic center and mean aerodynamic chord
    plt.scatter(MAC_vec[0], MAC_vec[1])
    plt.plot([MAC_LE_FS, MAC_TE_FS], [MAC_vec[1], MAC_vec[1]])

    for s in range(num_sections - 1):

        chord = section_Geo[s][0]
        FS_LE = planform_Geo[s][0]
        FS_TE = FS_LE + chord

        BL = planform_Geo[s][1]

        chord_next = section_Geo[s + 1][0]
        FS_LE_next = planform_Geo[s + 1][0]
        FS_TE_next = FS_LE_next + chord_next

        BL_next = planform_Geo[s + 1][1]

        LE_y = [FS_LE, FS_LE_next]
        BL_x = [BL, BL_next]
        neg_BL_x = [-BL, -BL_next]

        TE_y = [FS_TE, FS_TE_next]

        chord_y = [FS_LE, FS_TE]
        chord_x = [BL, BL]
        neg_chord_x = [-BL, -BL]

        #calculate hinge locations
        hinge_x = section_Geo[s][2]
        hinge_x_next = section_Geo[s + 1][2]
        hinge1_x_c = hinge_x * chord  # chordwise size of the hinge
        hinge2_x_c = hinge_x_next * chord_next  # chordwise size of next section
        hinge1_x = hinge1_x_c + FS_LE
        hinge2_x = hinge2_x_c + FS_LE_next
        hinge_index = planform_Geo[s][6]

        #color finder
        hinge_color = color_mat[int(hinge_index)]

        if hinge_index > -1:
            plt.plot([hinge1_x,hinge2_x],[BL, BL_next], c=hinge_color, linewidth=0.5)
            plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [BL, BL, BL_next, BL_next], c=hinge_color, alpha=0.4)
            if mirror == 1:
                plt.plot([hinge1_x, hinge2_x], [-BL, -BL_next], c=hinge_color, linewidth=0.5)
                plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [-BL, -BL, -BL_next, -BL_next], c=hinge_color, alpha=0.4)

        plt.plot(LE_y, BL_x, 'blue', linewidth=1)
        plt.plot(TE_y, BL_x, 'red', linewidth=1)
        plt.plot(chord_y, chord_x, 'black', linewidth=1)
        if mirror == 1:
            plt.plot(LE_y, neg_BL_x, 'blue', linewidth=1)
            plt.plot(TE_y, neg_BL_x, 'red', linewidth=1)
            plt.plot(chord_y, neg_chord_x, 'black', linewidth=1)

        if s == (num_sections - 2):
            chord_y = [FS_LE_next, FS_TE_next]
            chord_x = [BL_next, BL_next]
            neg_chord_x = [-BL_next, -BL_next]
            plt.plot(chord_y, chord_x, 'black', linewidth=1)
            if mirror == 1:
                plt.plot(chord_y,neg_chord_x, 'black', linewidth=1)

    plt.axis('equal')
    plt.grid()

    fig = plt.gcf()
    plt.xlabel('FS (X)')
    plt.ylabel('BL (Y)')

    return fig

#generate a matrix of colors for when a list of evenly spaced colors is needed
def gen_color_mat(min_color, max_color, num):

    if num <= 1:
        color_diff_0 = 0
        color_diff_1 = 0
        color_diff_2 = 0
    else:
        color_diff_0 = (max_color[0] - min_color[0]) / (num - 1)
        color_diff_1 = (max_color[1] - min_color[1]) / (num - 1)
        color_diff_2 = (max_color[2] - min_color[2]) / (num - 1)


    color_mat = [(0,0,0)]
    for i in range(num):
        c_offset = (i * color_diff_0, i * color_diff_1,
                          i * color_diff_2)
        hinge_color = (min_color[0] + c_offset[0], min_color[1] + c_offset[1],
                       min_color[2] + c_offset[2])
        color_mat.append(hinge_color)

    color_mat = color_mat[1::]

    return color_mat


def Plot_Planform_All(USER, save):

    plt.figure(figsize=(20, 20))

    body_list_path = gen_filename(USER, 'active_body_list', save)

    Body_List = process_string_csv(body_list_path, 0, 0)
    xlim = [0, 0]
    ylim = [0, 0]
    for body_i in range(len(Body_List)):

        Body = Body_List[body_i][0]

        body_json_filename = gen_filename(USER, f'{Body}_Geo.json', save)

        aero_Body = read_Aero_Body_from_json(body_json_filename)

        Body = aero_Body.name
        orientation = aero_Body.orientation
        body_mirror = aero_Body.mirror

        planform_Geo = aero_Body.planform_Geo
        section_Geo = aero_Body.section_Geo

        if orientation == 0:

            num_sections = len(section_Geo)
            num_hinge = max(planform_Geo[:, 6]) + 1
            if num_hinge == 0:
                num_hinge = 1

            color_mat = gen_color_mat((0.1, 0.1, 1), (1, 0.1, 0.25), int(num_hinge))

            MAC_vec = Calculate_Surface(aero_Body)
            MAC_vec = MAC_vec[0]

            MAC_LE_FS = (MAC_vec[2] * (-MAC_vec[3])) + MAC_vec[0]
            MAC_TE_FS = (MAC_vec[2] * (1 - MAC_vec[3])) + MAC_vec[0]

            #draw aerodynamic center and mean aerodynamic chord
            plt.scatter(MAC_vec[0], MAC_vec[1])
            plt.plot([MAC_LE_FS, MAC_TE_FS], [MAC_vec[1], MAC_vec[1]])

            for s in range(num_sections - 1):

                chord = section_Geo[s][0]
                FS_LE = planform_Geo[s][0]
                FS_TE = FS_LE + chord

                BL = planform_Geo[s][1]

                chord_next = section_Geo[s+1][0]
                FS_LE_next = planform_Geo[s+1][0]
                FS_TE_next = FS_LE_next + chord_next

                if FS_LE < xlim[0]:
                    xlim[0] = FS_LE

                if FS_TE_next > xlim[1]:
                    xlim[1] = FS_TE_next

                BL_next = planform_Geo[s+1][1]

                if BL_next > ylim[1]:
                    ylim = [-BL_next,BL_next]
                elif BL > ylim[1]:
                    ylim = [-BL,BL]

                LE_y = [FS_LE, FS_LE_next]
                BL_x = [BL, BL_next]
                neg_BL_x = [-BL, -BL_next]

                TE_y = [FS_TE, FS_TE_next]

                chord_y = [FS_LE, FS_TE]
                chord_x = [BL, BL]
                neg_chord_x = [-BL, -BL]

                # calculate hinge locations
                hinge_x = section_Geo[s][2]
                hinge_x_next = section_Geo[s + 1][2]
                hinge1_x_c = hinge_x * chord  # chordwise size of the hinge
                hinge2_x_c = hinge_x_next * chord_next  # chordwise size of next section
                hinge1_x = hinge1_x_c + FS_LE
                hinge2_x = hinge2_x_c + FS_LE_next
                hinge_index = planform_Geo[s][6]

                # color finder
                hinge_color = color_mat[int(hinge_index)]

                if hinge_index > -1:
                    plt.plot([hinge1_x, hinge2_x], [BL, BL_next], c=hinge_color, linewidth=0.5)
                    plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [BL,BL,BL_next,BL_next], c=hinge_color, alpha=0.4)
                    if body_mirror == 1:
                        plt.plot([hinge1_x, hinge2_x], [-BL, -BL_next], c=hinge_color, linewidth=0.5)
                        plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [-BL, -BL, -BL_next, -BL_next], c=hinge_color, alpha=0.4)

                plt.plot(LE_y, BL_x, 'blue', linewidth=1)
                plt.plot(TE_y, BL_x, 'red', linewidth=1)
                plt.plot(chord_y, chord_x, 'black', linewidth=1)
                if body_mirror == 1:
                    plt.plot(LE_y, neg_BL_x, 'blue', linewidth=1)
                    plt.plot(TE_y, neg_BL_x, 'red', linewidth=1)
                    plt.plot(chord_y, neg_chord_x, 'black', linewidth=1)

                if s == (num_sections - 2):
                    chord_y = [FS_LE_next, FS_TE_next]
                    chord_x = [BL_next, BL_next]
                    neg_chord_x = [-BL_next, -BL_next]
                    plt.plot(chord_y, chord_x, 'black', linewidth=1)
                    if body_mirror == 1:
                        plt.plot(chord_y,neg_chord_x, 'black', linewidth=1)

        elif orientation == 1:

            num_sections = len(section_Geo)

            for s in range(num_sections - 1):

                chord = section_Geo[s][0]
                FS_LE = planform_Geo[s][0]
                FS_TE = FS_LE + chord

                BL = planform_Geo[s][2]

                chord_next = section_Geo[s + 1][0]
                FS_LE_next = planform_Geo[s + 1][0]
                FS_TE_next = FS_LE_next + chord_next

                if FS_LE < xlim[0]:
                    xlim[0] = FS_LE

                if FS_TE_next > xlim[1]:
                    xlim[1] = FS_TE_next

                BL_next = planform_Geo[s + 1][2]

                if BL_next > ylim[1]:
                    ylim = [-BL_next,BL_next]
                elif BL > ylim[1]:
                    ylim = [-BL,BL]

                LE_y = [FS_LE, FS_LE_next]
                BL_x = [BL, BL_next]
                neg_BL_x = [-BL, -BL_next]

                TE_y = [FS_TE, FS_TE_next]

                chord_y = [FS_LE, FS_TE]
                chord_x = [BL, BL]
                neg_chord_x = [-BL, -BL]

                plt.plot(LE_y, BL_x, 'purple')
                plt.plot(TE_y, BL_x, 'yellow')
                plt.plot(chord_y, chord_x, 'black')

                if s == (num_sections - 2):
                    chord_y = [FS_LE_next, FS_TE_next]
                    chord_x = [BL_next, BL_next]
                    neg_chord_x = [-BL_next, -BL_next]
                    plt.plot(chord_y, chord_x, 'black')

    xlim = [xlim[0]-50, xlim[1]+50]
    ylim = [ylim[0] - 50, ylim[1] + 50]

    plt.axis('equal')
    plt.grid()
    plt.xlim(xlim)
    plt.ylim(ylim)

    fig = plt.gcf()

    plt.title('Top View')
    plt.xlabel('FS (X)')
    plt.ylabel('BL (Y)')

    return fig

def Plot_Planform_Allverts(USER, save):


    plt.figure(figsize=(20, 20))

    body_list_path = gen_filename(USER, 'active_body_list', save)

    Body_List = process_string_csv(body_list_path, 0, 0)
    xlim = [0, 0]
    ylim = [0, 0]
    for body_i in range(len(Body_List)):

        Body = Body_List[body_i][0]

        body_json_filename = gen_filename(USER, f'{Body}_Geo.json', save)

        aero_Body = read_Aero_Body_from_json(body_json_filename)

        Body = aero_Body.name
        orientation = aero_Body.orientation
        body_mirror = aero_Body.mirror

        planform_Geo = aero_Body.planform_Geo
        section_Geo = aero_Body.section_Geo

        if orientation == 1:

            num_sections = len(section_Geo)
            num_hinge = max(planform_Geo[:, 6]) + 1
            if num_hinge == 0:
                num_hinge = 1

            color_mat = gen_color_mat((0.1, 0.1, 1), (1, 0.1, 0.25), int(num_hinge))

            MAC_vec = Calculate_Surface(aero_Body)
            MAC_vec = MAC_vec[0]

            MAC_LE_FS = (MAC_vec[2] * (-MAC_vec[3])) + MAC_vec[0]
            MAC_TE_FS = (MAC_vec[2] * (1 - MAC_vec[3])) + MAC_vec[0]

            #draw aerodynamic center and mean aerodynamic chord
            plt.scatter(MAC_vec[0], MAC_vec[1])
            plt.plot([MAC_LE_FS, MAC_TE_FS], [MAC_vec[1], MAC_vec[1]])

            for s in range(num_sections - 1):

                chord = section_Geo[s][0]
                FS_LE = planform_Geo[s][0]
                FS_TE = FS_LE + chord

                BL = planform_Geo[s][1]

                chord_next = section_Geo[s+1][0]
                FS_LE_next = planform_Geo[s+1][0]
                FS_TE_next = FS_LE_next + chord_next

                if FS_LE < xlim[0]:
                    xlim[0] = FS_LE

                if FS_TE_next > xlim[1]:
                    xlim[1] = FS_TE_next

                BL_next = planform_Geo[s+1][1]

                if BL_next > ylim[1]:
                    ylim = [-BL_next,BL_next]
                elif BL > ylim[1]:
                    ylim = [-BL,BL]

                LE_y = [FS_LE, FS_LE_next]
                BL_x = [BL, BL_next]
                neg_BL_x = [-BL, -BL_next]

                TE_y = [FS_TE, FS_TE_next]

                chord_y = [FS_LE, FS_TE]
                chord_x = [BL, BL]
                neg_chord_x = [-BL, -BL]

                # calculate hinge locations
                hinge_x = section_Geo[s][2]
                hinge_x_next = section_Geo[s + 1][2]
                hinge1_x_c = hinge_x * chord  # chordwise size of the hinge
                hinge2_x_c = hinge_x_next * chord_next  # chordwise size of next section
                hinge1_x = hinge1_x_c + FS_LE
                hinge2_x = hinge2_x_c + FS_LE_next
                hinge_index = planform_Geo[s][6]

                # color finder
                hinge_color = color_mat[int(hinge_index)]

                if hinge_index > -1:
                    plt.plot([hinge1_x, hinge2_x], [BL, BL_next], c=hinge_color, linewidth=0.5)
                    plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [BL, BL, BL_next, BL_next], c=hinge_color, alpha=0.4)
                    if body_mirror == 1:
                        plt.plot([hinge1_x, hinge2_x], [-BL, -BL_next], c=hinge_color, linewidth=0.5)
                        plt.fill([FS_TE, hinge1_x, hinge2_x, FS_TE_next], [-BL, -BL, -BL_next, -BL_next], c=hinge_color, alpha=0.4)

                plt.plot(LE_y, BL_x, 'purple', linewidth=1)
                plt.plot(TE_y, BL_x, 'yellow', linewidth=1)
                plt.plot(chord_y, chord_x, 'black', linewidth=1)
                if body_mirror == 1:
                    plt.plot(LE_y, neg_BL_x, 'purple', linewidth=1)
                    plt.plot(TE_y, neg_BL_x, 'yellow', linewidth=1)
                    plt.plot(chord_y, neg_chord_x, 'black', linewidth=1)

                if s == (num_sections - 2):
                    chord_y = [FS_LE_next, FS_TE_next]
                    chord_x = [BL_next, BL_next]
                    neg_chord_x = [-BL_next, -BL_next]
                    plt.plot(chord_y, chord_x, 'black', linewidth=1)
                    if body_mirror == 1:
                        plt.plot(chord_y,neg_chord_x, 'black', linewidth=1)

        elif orientation == 0:

            num_sections = len(section_Geo)

            for s in range(num_sections - 1):

                chord = section_Geo[s][0]
                FS_LE = planform_Geo[s][0]
                FS_TE = FS_LE + chord

                BL = planform_Geo[s][2]

                chord_next = section_Geo[s + 1][0]
                FS_LE_next = planform_Geo[s + 1][0]
                FS_TE_next = FS_LE_next + chord_next

                if FS_LE < xlim[0]:
                    xlim[0] = FS_LE

                if FS_TE_next > xlim[1]:
                    xlim[1] = FS_TE_next

                BL_next = planform_Geo[s + 1][2]

                if BL_next > ylim[1]:
                    ylim = [-BL_next,BL_next]
                elif BL > ylim[1]:
                    ylim = [-BL,BL]

                LE_y = [FS_LE, FS_LE_next]
                BL_x = [BL, BL_next]
                neg_BL_x = [-BL, -BL_next]

                TE_y = [FS_TE, FS_TE_next]

                chord_y = [FS_LE, FS_TE]
                chord_x = [BL, BL]
                neg_chord_x = [-BL, -BL]

                plt.plot(LE_y, BL_x, 'blue')
                plt.plot(TE_y, BL_x, 'red')
                plt.plot(chord_y, chord_x, 'black')

                if s == (num_sections - 2):
                    chord_y = [FS_LE_next, FS_TE_next]
                    chord_x = [BL_next, BL_next]
                    neg_chord_x = [-BL_next, -BL_next]
                    plt.plot(chord_y, chord_x, 'black')

    plt.axis('equal')
    plt.grid()

    xlim = [xlim[0] - 50, xlim[1] + 50]
    ylim = [ylim[0] - 50, ylim[1] + 50]

    plt.xlim(xlim)
    plt.ylim(ylim)

    fig = plt.gcf()

    plt.title('Side View')
    plt.xlabel('FS (X)')
    plt.ylabel('WL (Z)')

    return fig
