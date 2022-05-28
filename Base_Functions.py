#Basic level functions that are generic and not specific to vortex master

import numpy as np
import matplotlib.pyplot as plt


#change from degrees to radians
def deg2rad(angle_Deg):

    angle_in_rad = angle_Deg * (np.pi / 180)

    return angle_in_rad

def rad2deg(angle_Rad):

    angle_in_deg = angle_Rad * (180 / np.pi)

    return  angle_in_deg

#rotate a 2-D array of points a specified amount about a point
def rotate_about_point(Point_List, rotation_angle, rotation_x, rotation_z):

    translated2hinge = Point_List[:,0] - rotation_x, Point_List[:,1] - rotation_z

    rotatedx = translated2hinge[0] * np.cos(deg2rad(-rotation_angle)) - translated2hinge[1] * np.sin(deg2rad(-rotation_angle))
    rotatedz = translated2hinge[1] * np.cos(deg2rad(-rotation_angle)) + translated2hinge[0] * np.sin(deg2rad(-rotation_angle))

    translated_back = rotatedx + rotation_x, rotatedz + rotation_z

    return translated_back

def rotate_surface_list(surface_list_x, surface_list_y, rotation_angle, rotation_x, rotation_y):

    cos = np.cos(deg2rad(-rotation_angle))
    sin = np.sin(deg2rad(-rotation_angle))

    for i in range(len(surface_list_x)):


        for n in range(len(surface_list_x[i])):

            x = surface_list_x[i][n]
            y = surface_list_y[i][n]

            temp_x = x - rotation_x
            temp_y = y - rotation_y

            temp_x_1 = (temp_x * cos) - (temp_y * sin)
            temp_y_1 = (temp_x * sin) + (temp_y * cos)

            new_x = temp_x_1 + rotation_x
            new_y = temp_y_1 + rotation_y

            surface_list_x[i][n] = new_x
            surface_list_y[i][n] = new_y

    return [surface_list_x, surface_list_y]

def circle_pts(center_x, center_z, radius):

    theta = np.linspace(0, 2 * np.pi, 100)

    x = radius * np.cos(theta) + center_x
    z = radius * np.sin(theta) + center_z

    circle_array = np.array([x, z])

    return circle_array

def arc_pts(center_x, center_y, radius, start_angle, end_angle, num_pts):

    theta = np.linspace(start_angle, end_angle, num_pts)

    x = radius * np.cos(theta) + center_x
    y = radius * np.sin(theta) + center_y


    return x, y

def rotate_scale_translate(pts, scale, rotation_angle, rotation_x, rotation_z, translate_x, translate_z):

    rotated_pts = rotate_about_point(pts, rotation_angle, rotation_x,rotation_z)

    scaled_pts = np.multiply(rotated_pts, scale)

    scaled_pts[0] = scaled_pts[0] + translate_x
    scaled_pts[1] = scaled_pts[1] + translate_z

    return scaled_pts

#check string for any commas or spaces
def str_check(str):

    space_check = any(i.isspace() for i in str)
    comma_checc = str.find(',')

    if comma_checc == -1:
        comma_check = 1
    else:
        comma_check = 0

    if space_check == True:
        space_check = 0
    else:
        space_check = 1

    comb_check = comma_check * space_check

    print(comb_check)

    return comb_check

def tri_angle_from_sides(c, b, a):

    C = np.arccos((a**2 + b**2 - c**2) / (2 * a * b))
    C = rad2deg(C)

    B = np.arccos((a**2 + c**2 - b**2) / (2 * a * c))
    B = rad2deg(B)

    A = np.arccos((b**2 + c**2 - a**2) / (2 * b * c))
    A = rad2deg(A)

    return [C, B, A]

def pt_pt_distance(pt1_x, pt1_y, pt2_x, pt2_y):

    diff_x = pt2_x- pt1_x
    diff_y = pt2_y- pt1_y

    dis = np.sqrt(diff_x ** 2 + diff_y ** 2).tolist()

    return dis