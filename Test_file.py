import time

from Data_Managment import *
from Data_Processing import *
from Plotting import *
from Base_Calculations import *

import matplotlib.pyplot as plt
import json
from numpy import ndarray

user = 'user2'
save_name = 'aircraft1'
body_name = 'Wing'

body_json_filename = gen_filename(user, f'{body_name}_Geo.json', save_name)

aero_Body = read_Aero_Body_from_json(body_json_filename)


new_airfoil = Airfoil_Section(aero_Body.section_Geo[0], aero_Body.airfoil_Geo[0][0])


#plt.plot(new_airfoil.upper_x,new_airfoil.upper_y)
#plt.plot(new_airfoil.lower_x,new_airfoil.lower_y)
#plt.scatter(new_airfoil.upper_x[new_airfoil.min_dis_index],new_airfoil.upper_y[new_airfoil.min_dis_index])
#plt.scatter(new_airfoil.upper_x[new_airfoil.min_dis_index2],new_airfoil.upper_y[new_airfoil.min_dis_index2])

#plt.scatter(new_airfoil.lower_x[new_airfoil.min_low_index],new_airfoil.lower_y[new_airfoil.min_low_index])
#plt.scatter(new_airfoil.lower_x[new_airfoil.min_low_index2],new_airfoil.lower_y[new_airfoil.min_low_index2])

circ = circle_pts(new_airfoil.hinge_x, new_airfoil.hinge_y, new_airfoil.min_hinge_rad)
plt.plot(circ[0], circ[1],'--', linewidth=1)

hinge_rotation = 25 #in degrees
[new_airfoil.hinge_surfaces_x, new_airfoil.hinge_surfaces_y] = rotate_surface_list(new_airfoil.hinge_surfaces_x, new_airfoil.hinge_surfaces_y,hinge_rotation,new_airfoil.hinge_x, new_airfoil.hinge_y)


for n in range(len(new_airfoil.bottom_surfaces_x)):

    x = new_airfoil.bottom_surfaces_x[n]
    y = new_airfoil.bottom_surfaces_y[n]

    plt.plot(x, y, linewidth=1, color='red')

for n in range(len(new_airfoil.upper_surfaces_x)):
    x = new_airfoil.upper_surfaces_x[n]
    y = new_airfoil.upper_surfaces_y[n]

    plt.plot(x, y, linewidth=1, color='blue')

for n in range(len(new_airfoil.hinge_surfaces_x)):
    x = new_airfoil.hinge_surfaces_x[n]
    y = new_airfoil.hinge_surfaces_y[n]

    plt.plot(x, y, linewidth=1, color='green')

#plt.plot(new_airfoil.hinge_lower_x, new_airfoil.hinge_lower_y, linewidth=1)

#plt.plot(new_airfoil.hinge_pts_x, new_airfoil.hinge_pts_y)

#plt.scatter(new_airfoil.hinge_x, new_airfoil.hinge_y)
#plt.scatter(new_airfoil.meet_1x, new_airfoil.meet_1y)
#plt.scatter(new_airfoil.meet_2x, new_airfoil.meet_2y)

plt.axis('equal')
plt.show()


plt.Circle((new_airfoil.hinge_x, new_airfoil.hinge_y), new_airfoil.min_hinge_rad)

#fig = Plot_Planform(aero_Body)
#plt.savefig(f'user_data\\{user}\\images\\body_planform.jpg')

fig = Plot_Planform_All(user, save_name)
#plt.savefig(f'user_data\\{user}\\images\\all_planform.jpg')

fig = Plot_Planform_Allverts(user, save_name)
#plt.savefig(f'user_data\\{user}\\images\\all_vert_planform.jpg')

Plotting_data = process_plot_body_object(aero_Body)

fig = Plot_Airfoil_Section(0,Plotting_data)

json_filename = 'user_data\\user2\\files\\aircraft1\\temp_save_data\\Canard_Struct.json'

#plt.show()











