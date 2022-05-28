#functions to interact with save and temp files
#also manage log files

from Data_Processing import *
from Base_Functions import *
import json


inertial_headers = ['item','mass','CG_FS','CG_BL','CG_WL']
balance_headers = ['FS', 'BL', 'WL', 'ac', 'Total Area']

planform_headers = ['LE_FS', 'BL', 'WL', 'Section Span', 'LE_Sweep', 'Dihedral', 'Hinge_index', 'Input Type']
section_headers = ['Chord', 'Twist', 'Hinge x/c', 'Hinge y/c', 'Hinge r/c']
airfoil_headers = ['Airfoils']

Body_headers = ['Body name', 'orientation', 'mirrored', 'num_hinge']
save_Body_headers = ['Save Info', 'Aircraft Name']

#object that holds all basic geometry of a body
class Aero_Body:
    def __init__(self, name, orientation, mirror, num_hinge, planform_Geo, section_Geo, airfoil_Geo):
        self.name = name
        self.orientation = orientation
        self.mirror = mirror
        self.num_hinge = num_hinge

        self.planform_Geo = planform_Geo
        self.section_Geo = section_Geo
        self.airfoil_Geo = airfoil_Geo

#defines a single cross section to calculate hinge size and other variables
class Airfoil_Section:
    def __init__(self, section_row, airfoil):
        self.section_row = section_row
        self.airfoil = str(airfoil)

        self.solve_geometry()
        self.solve_new_profiles()

    #gets the list for upper and lower surfaces and calculates the minimum hinge radius based on hinge location
    def solve_geometry(self):
        [airfoil_database_strings, airfoil_database_x, airfoil_database_y] = get_airfoil_database2()

        airfoil_index = airfoil_database_strings.index(self.airfoil)

        airfoil_x = airfoil_database_x[airfoil_index]
        airfoil_y = airfoil_database_y[airfoil_index]


        #seperate upper and lower sections of airfoil
        upper_x = []
        upper_y = []
        lower_x = [0] #give initial origin point
        lower_y = [0] #give initial origin point
        on_up = 1
        for i in range(len(airfoil_x)):

            x = airfoil_x[i]
            y = airfoil_y[i]

            if on_up == 1:
                upper_x.append(x)
                upper_y.append(y)

            elif on_up == 0:
                lower_x.append(x)
                lower_y.append(y)

            if i > 0 and x == 0 and y == 0:
                on_up = 0

        self.upper_x = upper_x
        self.upper_y = upper_y
        self.lower_x = lower_x
        self.lower_y = lower_y


        #find min radius of hinge circle to begin design of hinge profile
        hinge_x = self.section_row[2]
        hinge_y = self.section_row[3]

        self.hinge_x = hinge_x
        self.hinge_y = hinge_y

        hinge_diff_x = self.upper_x - hinge_x
        hinge_diff_y = self.upper_y - hinge_y

        hinge_dis = np.sqrt(hinge_diff_x**2 + hinge_diff_y**2).tolist() #list of distance to each point on upper surface

        min_dis_index = hinge_dis.index(min(hinge_dis))

        self.min_dis_index = min_dis_index

        #find second lowest distance
        hinge_dis2 = hinge_dis
        hinge_dis2[self.min_dis_index] = 1
        min_dis_index2 = hinge_dis2.index(min(hinge_dis2))

        self.min_dis_index2 = min_dis_index2

        #solve for interpolated closest point between 2 lowest distance given points
        c = np.sqrt((upper_x[min_dis_index] - upper_x[min_dis_index2])**2 + (upper_y[min_dis_index] - upper_y[min_dis_index2])**2)

        b = np.sqrt((hinge_x - upper_x[min_dis_index])**2 + (hinge_y - upper_y[min_dis_index])**2)
        a = np.sqrt((hinge_x - upper_x[min_dis_index2]) ** 2 + (hinge_y - upper_y[min_dis_index2]) ** 2)

        [C, B, A] = tri_angle_from_sides(c, b, a)

        C1 = 90 - A

        l1_c = b * np.cos(deg2rad(A)) / c

        #meeting point of hinge against upper surface
        self.meet_1x = upper_x[min_dis_index] - (l1_c * (upper_x[min_dis_index] - upper_x[min_dis_index2]))
        self.meet_1y = upper_y[min_dis_index] - (l1_c * (upper_y[min_dis_index] - upper_y[min_dis_index2]))

        self.min_hinge_rad = b * np.cos(deg2rad(C1))

        low_base_check_x = []
        low_base_check_y = []
        for n in range(len(lower_x)):

            x = lower_x[n]
            y = lower_y[n]

            if x > self.hinge_x:
                break
            else:
                low_base_check_x.append(x)
                low_base_check_y.append(y)

        low_hinge_diff_x = low_base_check_x - hinge_x
        low_hinge_diff_y = low_base_check_y - hinge_y

        low_hinge_dis = np.sqrt(low_hinge_diff_x**2 + low_hinge_diff_y**2).tolist()

        low_dis_from_radius = abs(low_hinge_dis - self.min_hinge_rad).tolist()

        self.min_low_index = low_dis_from_radius.index(min(low_dis_from_radius))

        low_dis_from_radius2 = low_dis_from_radius
        low_dis_from_radius2[self.min_low_index] = 1
        min_low_index2 = low_dis_from_radius2.index(min(low_dis_from_radius2))

        self.min_low_index2 = min_low_index2

        dis1 = low_hinge_dis[self.min_low_index]
        dis2 = low_hinge_dis[self.min_low_index2]

        perc = (self.min_hinge_rad - dis1) / (dis2 - dis1)

        self.meet_2x = lower_x[self.min_low_index] - (perc * (lower_x[self.min_low_index] - lower_x[self.min_low_index2]))
        self.meet_2y = lower_y[self.min_low_index] - (perc * (lower_y[self.min_low_index] - lower_y[self.min_low_index2]))

    #solve for new surfaces accounting for hinge
    #NEED TO ACCOUNT FOR GAP IN HINGE AND SMOOTH TRANSITION TO HINGE ON BOTTOM SURFACE
    def solve_new_profiles(self):

        base_upper_x = [self.meet_1x]
        hinge_upper_x = []
        base_upper_y = [self.meet_1y]
        hinge_upper_y = []

        on_base = 0
        for i in range(len(self.upper_x)):

            x = self.upper_x[i]
            y = self.upper_y[i]

            if i > 0 and x < self.meet_1x:

                on_base = 1
                hinge_upper_x.append(self.meet_1x)
                hinge_upper_y.append(self.meet_1y)

            if on_base == 0:

                hinge_upper_x.append(x)
                hinge_upper_y.append(y)

            elif on_base == 1:

                base_upper_x.append(x)
                base_upper_y.append(y)

        self.base_upper_x = base_upper_x
        self.hinge_upper_x = hinge_upper_x
        self.base_upper_y = base_upper_y
        self.hinge_upper_y = hinge_upper_y

        base_lower_x = []
        hinge_lower_x = [self.meet_2x]
        base_lower_y = []
        hinge_lower_y = [self.meet_2y]

        on_base = 1
        for n in range(len(self.lower_x)):

            x = self.lower_x[n]
            y = self.lower_y[n]

            if n > 0 and x > self.meet_2x:
                on_base = 0

                base_lower_x.append(self.meet_2x)
                base_lower_y.append(self.meet_2y)

            if on_base == 1:

                base_lower_x.append(x)
                base_lower_y.append(y)

            elif on_base == 0:

                hinge_lower_x.append(x)
                hinge_lower_y.append(y)

        self.base_lower_x = base_lower_x
        self.hinge_lower_x = hinge_lower_x
        self.base_lower_y = base_lower_y
        self.hinge_lower_y = hinge_lower_y

        meet1_2 = pt_pt_distance(self.meet_1x, self.meet_1y, self.meet_2x, self.meet_2y)

        [b, a, MEETh] = tri_angle_from_sides(self.min_hinge_rad, self.min_hinge_rad, meet1_2)

        abs_deg_h_meet1 = np.arctan((self.meet_1y - self.hinge_y) / (self.meet_1x - self.hinge_x))

        num_hinge_pts = 50

        [hinge_pts_x, hinge_pts_y] = arc_pts(self.hinge_x, self.hinge_y, self.min_hinge_rad, abs_deg_h_meet1,
                                             abs_deg_h_meet1 + deg2rad(MEETh), num_hinge_pts)

        self.hinge_pts_x = hinge_pts_x.tolist()
        self.hinge_pts_y = hinge_pts_y.tolist()

        self.hinge_pts_x.reverse()
        self.hinge_pts_y.reverse()


        self.bottom_surfaces_x = [self.base_lower_x, hinge_pts_x]
        self.bottom_surfaces_y = [self.base_lower_y, hinge_pts_y]

        self.upper_surfaces_x = [self.base_upper_x]
        self.upper_surfaces_y = [self.base_upper_y]

        self.hinge_surfaces_x = [self.hinge_upper_x, self.hinge_pts_x, self.hinge_lower_x]
        self.hinge_surfaces_y = [self.hinge_upper_y, self.hinge_pts_y, self.hinge_lower_y]


class Struct_Body:
    def __init__(self, name, name_list, mass_Geo, mass_Size, struct_Geo):
        self.name = name
        self.name_list = name_list
        self.mass_Geo = mass_Geo
        self.mass_Size = mass_Size
        self.struct_Geo = struct_Geo

        self.solve_inertial()

    #solves the average inertial properties of the body
    def solve_inertial(self):
        self.total_mass = float(sum(np.array(self.mass_Geo)[:,0]))

        sum_CG_FS = 0
        sum_CG_BL = 0
        sum_CG_WL = 0

        sum_Ixx = 0
        sum_Iyy = 0
        sum_Izz = 0
        sum_Iyz = 0
        sum_Ixz = 0
        sum_Ixy = 0

        num_items = len(self.name_list)

        for n in range(num_items):
            row = self.mass_Geo[n]

            mass = row[0]
            item_FS = row[1]
            item_BL = row[2]
            item_WL = row[3]

            sum_CG_FS += mass * item_FS
            sum_CG_BL += mass * item_BL
            sum_CG_WL += mass * item_WL

        CG_FS = sum_CG_FS / self.total_mass
        CG_BL = sum_CG_BL / self.total_mass
        CG_WL = sum_CG_WL / self.total_mass

        for i in range(num_items):
            row = self.mass_Geo[i]

            mass = row[0]
            item_FS = row[1]
            item_BL = row[2]
            item_WL = row[3]

            sum_Ixx += mass * ((item_BL - CG_BL) ** 2 + (item_WL - CG_WL) ** 2)
            sum_Iyy += mass * ((item_FS - CG_FS) ** 2 + (item_WL - CG_WL) ** 2)
            sum_Izz += mass * ((item_FS - CG_FS) ** 2 + (item_BL - CG_BL) ** 2)
            sum_Iyz += mass * ((item_BL - CG_BL) * (item_WL - CG_WL))
            sum_Ixz += mass * ((item_FS - CG_FS) * (item_WL - CG_WL))
            sum_Ixy += mass * ((item_FS - CG_FS) * (item_BL - CG_BL))

        self.Ixx = sum_Ixx
        self.Iyy = sum_Iyy
        self.Izz = sum_Izz
        self.Iyz = sum_Iyz
        self.Ixz = sum_Ixz
        self.Ixy = sum_Ixy

#find filenames of body files based on user, save, and body
def get_filenames(user, save_name, body_name):

    Planform_filename = 'user_data\{0}\\files\{1}\\temp_save_data\Geometry_Files\{2}_Planform_Geometry.csv'.format(user, save_name, body_name)
    Section_filename = 'user_data\{0}\\files\{1}\\temp_save_data\Geometry_Files\{2}_Section_Geometry.csv'.format(user, save_name, body_name)
    Airfoil_filename = 'user_data\{0}\\files\{1}\\temp_save_data\Geometry_Files\{2}_Airfoil_Geometry.csv'.format(user, save_name, body_name)

    Save_list_filename = 'user_data\{0}\\active_file_list.csv'.format(user, save_name, body_name)

    return [Planform_filename, Section_filename, Airfoil_filename, Save_list_filename]

#write planform, section, and airfoil to a .json file for js code to pick up
def write_body_to_json(json_filename, planform_Geo, section_Geo, airfoil_Geo):

    planform_Geo = planform_Geo.tolist()
    section_Geo = section_Geo.tolist()
    airfoil_Geo = airfoil_Geo.tolist()

    body_data = Body_Data(planform_Geo, section_Geo, airfoil_Geo)

    with open(json_filename, 'w') as outfile:
        json.dump(body_data.__dict__, outfile, indent=4)

#write planform, section, and airfoil to a .json file for js code to pick up
def write_Aero_Body_to_json(json_filename, name, orientation, mirror, num_hinge, planform_Geo, section_Geo, airfoil_Geo):

    name = str(name)
    orientation = int(orientation)
    mirror = int(mirror)
    num_hinge = int(num_hinge)

    planform_Geo = planform_Geo.tolist()
    section_Geo = section_Geo.tolist()
    airfoil_Geo = airfoil_Geo.tolist()

    aero_Body = Aero_Body(name, orientation, mirror, num_hinge, planform_Geo, section_Geo, airfoil_Geo)

    with open(json_filename, 'w') as outfile:
        json.dump(aero_Body.__dict__, outfile, indent=4)


#saves an object data type into a .json file
def update_json_file(filename, object_Type):

    with open(filename, 'w') as outfile:
        json.dump(object_type.__dict__, outfile, indent=4)

#read the body data from json file to use in python backend
def read_Aero_Body_from_json(body_json_filename):

    with open(body_json_filename) as json_file:
        data = json.load(json_file)

    name = data['name']
    orientation = data['orientation']
    mirror = data['mirror']
    num_hinge = data['num_hinge']

    planform_Geo = np.array(data['planform_Geo'])
    section_Geo = np.array(data['section_Geo'])
    airfoil_Geo = np.array(data['airfoil_Geo'])

    aero_Body = Aero_Body(name, orientation, mirror, num_hinge, planform_Geo, section_Geo, airfoil_Geo)

    return aero_Body

#read the body data from json file to use in python backend
def read_Struct_Body_from_json(struct_json_filename):

    with open(struct_json_filename) as json_file:
        data = json.load(json_file)

    name = data['name']

    name_list = data['name_list']
    mass_Geo = data['mass_Geo']
    mass_Size = data['mass_Size']
    struct_Geo = data['struct_Geo']

    struct_Body = Struct_Body(name, name_list, mass_Geo, mass_Size, struct_Geo)

    return struct_Body

#read the body data from json file to use in python backend
def read_body_from_json(body_json_filename):

    with open(body_json_filename) as json_file:
        data = json.load(json_file)

    planform_Geo = np.array(data['planform_Geo'])
    section_Geo = np.array(data['section_Geo'])
    airfoil_Geo = np.array(data['airfoil_Geo'])

    body_data = Body_Data(planform_Geo,section_Geo,airfoil_Geo)

    return [planform_Geo, section_Geo, airfoil_Geo, body_data]

#generates a filename string for any file needed for data management
#designed to make file structure more managable and easier to change
def gen_filename(USER, des_file, save):

    if USER == 'SHARED':
        user_dir = 'shared_data\\'
    else:
        user_dir = f'user_data\{USER}\\'

    if des_file == 'user_preferences':
        type_dir = ''
    elif des_file == 'active_file_list':
        type_dir = ''
    elif 'database' in des_file:
        type_dir = 'databases\\'
    elif 'log' in des_file:
        type_dir = 'logs\\'
    elif 'Geometry' in des_file:
        type_dir = f'files\{save}\\temp_save_data\geometry_files\\'
    elif 'Geo' in des_file:
        type_dir = f'files\{save}\\temp_save_data\\'
    elif 'Properties' in des_file:
        type_dir = f'files\{save}\\temp_save_data\\'
    elif des_file == 'active_body_list':
        type_dir = f'files\{save}\\temp_save_data\\'

    if '.jpg' in des_file:
        file_type = ''
    if '.json' in des_file:
        file_type = ''
    else:
        file_type = '.csv'

    filename = f'{user_dir}{type_dir}{des_file}{file_type}'

    return filename

#fuction to get all airfoil database info
#in later version, merge user data with shared data
def get_airfoil_database():

    shared_airfoil_database_path = gen_filename('SHARED', 'airfoil_database', 0)

    Airfoil_Database = process_csv(shared_airfoil_database_path, 3, 0)
    Airfoil_Database_strings = process_string_csv(shared_airfoil_database_path, 0, 0)
    Airfoil_Database_strings = Airfoil_Database_strings[0, :]
    Airfoil_Database_strings = make_proper_list(Airfoil_Database_strings)

    return [Airfoil_Database, Airfoil_Database_strings]

def get_airfoil_database2():

    shared_airfoil_database_path = gen_filename('SHARED', 'airfoil_database.json', 0)

    with open(shared_airfoil_database_path) as json_file:
        data = json.load(json_file)

    airfoil_database_strings = data['name_list']
    airfoil_database_x = data['x_list']
    airfoil_database_y = data['y_list']

    return [airfoil_database_strings, airfoil_database_x, airfoil_database_y]


#rewrite airfoil database with new airfoil added to end
def add_airfoil_to_list(database_filename, new_airfoil_mat):

    main_airfoil_array = process_string_csv(database_filename, 0, 0)

    len_main = len(main_airfoil_array)
    len_new = len(new_airfoil_mat)

    if len_main > len_new:

        extra = ['', '']

        for i in range(len_main - len_new):

            extra = np.vstack((extra, ['1','0']))

        extra = extra[1::]

        equalized_new_array = np.vstack((new_airfoil_mat, extra))

        new_full = np.hstack((main_airfoil_array, equalized_new_array))

    elif len_new >= len_main:

        new_airfoil_mat = np.vstack((new_airfoil_mat, ['1','0']))

        len_new = len_new + 1

        num_col = len(main_airfoil_array[0])

        extra = np.zeros(num_col)

        full_row = []

        x = num_col/2
        x = np.int(x)

        for i in range(x):

            full_row = np.hstack((full_row, ['1','0']))

        for i in range(len_new - len_main):

            extra = np.vstack((extra, full_row))

        extra = extra[1::]

        equalized_main = np.vstack((main_airfoil_array, extra))

        new_full = np.hstack((equalized_main, new_airfoil_mat))

    num_col = len(new_full[0])

    initial_zeros = []
    for n in range(num_col):
        initial_zeros.append('0')

    new_full = np.vstack((initial_zeros, new_full))

    file = open(database_filename, 'w+', newline ='')
    with file:
        write = csv.writer(file)
        write.writerows(new_full)

#rewrite the geometry files with the values from the user input matrix
def update_body_files(values, filename, callback_mat, type):

    current_array = process_string_csv(filename, 0, 0)

    if type == 'planform':
        headers = planform_headers
        input_type_list = current_array[:, 7]
    elif type == 'section':
        headers = section_headers
    elif type == 'airfoil':
        headers = airfoil_headers

    new_mat = [headers]

    for r in range(len(callback_mat)):

        new_row = []

        for c in range(len(callback_mat[0])):
            callback = callback_mat[r][c]

            cell_val = values[callback]

            if cell_val == 'Position':
                cell_val = 0
            elif cell_val == 'Angles':
                cell_val = 1

            new_row.append(cell_val)

        new_mat.append(new_row)

    file = open(filename, 'w+', newline='')
    with file:
        write = csv.writer(file)
        write.writerows(new_mat)

#shouldnt be needed in js version
def update_file_from_GUItable(values, filename, callback_mat, headers):

    new_mat = [headers]

    for r in range(len(callback_mat)):

        new_row = []

        for c in range(len(callback_mat[0])):
            callback = callback_mat[r][c]

            cell_val = values[callback]

            new_row.append(cell_val)

        new_mat.append(new_row)

    file = open(filename, 'w+', newline='')
    with file:
        write = csv.writer(file)
        write.writerows(new_mat)

#rewrite the geometry files with the values from the user input matrix
def update_body_files_with_mat(new_mat, filename, type):

    current_array = process_string_csv(filename, 0, 0)

    if type == 'planform':
        headers = planform_headers
        input_type_list = current_array[:, 7]
    elif type == 'section':
        headers = section_headers
    elif type == 'airfoil':
        headers = airfoil_headers

    new_full_mat = [headers]

    new_full_mat = np.vstack((new_full_mat,new_mat))

    file = open(filename, 'w+', newline='')
    with file:
        write = csv.writer(file)
        write.writerows(new_full_mat)

def update_file_with_mat(new_mat, filename, headers):
    new_full_mat = [headers]

    new_full_mat = np.vstack((new_full_mat, new_mat))

    file = open(filename, 'w+', newline='')
    with file:
        write = csv.writer(file)
        write.writerows(new_full_mat)

#add a section in a specified spot in geometry arrays
def add_section_to_body(aero_Body, user, save_name, section_index, add_at_percentage):

    body_name = aero_Body.name
    planform_Geo = aero_Body.planform_Geo
    section_Geo = aero_Body.section_Geo
    airfoil_Geo = aero_Body.airfoil_Geo

    if int(section_index) == len(section_Geo) - 1:
        print('Cannot add section after last section')
    else:
        add_at_percentage = np.float(add_at_percentage)

        prev_Planform_mat = planform_Geo[:section_index + 1]
        next_Planform_mat = planform_Geo[section_index + 1::]

        prev_Planform_row = prev_Planform_mat[-1]
        next_Planform_row = next_Planform_mat[0]
        new_Planform_row = [0]

        new_Planform_mat = [[0,0,0,0,0,0,0,0]]
        for n in range(len(prev_Planform_row)):

            prev_val = np.float(prev_Planform_row[n])
            next_val = np.float(next_Planform_row[n])

            #fs,bl,wl averaging
            if n <= 2:
                new_val = (prev_val + next_val) * add_at_percentage

            elif n == 3:
                new_val = prev_Planform_row[n] * (1 - add_at_percentage)

                new_prev_val = prev_val - new_val
                prev_Planform_row[n] = new_prev_val

            elif n >= 4 and n <= 6:

                new_val = prev_Planform_row[n]

            elif n == 7:

                new_val = 1

            new_Planform_row.append(new_val)

        new_Planform_row = new_Planform_row[1::]

        prev_Planform_mat[-1] = prev_Planform_row
        next_Planform_mat[0] = next_Planform_row

        new_Planform_mat = np.vstack((new_Planform_mat,prev_Planform_mat))
        new_Planform_mat = np.vstack((new_Planform_mat,new_Planform_row))
        new_Planform_mat = np.vstack((new_Planform_mat,next_Planform_mat))

        new_Planform_mat = new_Planform_mat[1::]

        ##################################

        prev_Section_mat = section_Geo[:section_index + 1]
        next_Section_mat = section_Geo[section_index + 1::]

        prev_Section_row = prev_Section_mat[-1]
        next_Section_row = next_Section_mat[0]
        new_Section_row = [0]

        new_Section_mat = [[0, 0, 0, 0, 0]]
        for n in range(len(prev_Section_row)):

            prev_val = np.float(prev_Section_row[n])
            next_val = np.float(next_Section_row[n])

            if n == 1:
                new_val = (prev_val + next_val) * add_at_percentage
            elif n == 0:
                new_val = prev_val + ((next_val - prev_val) * add_at_percentage)
            else:
                new_val = (prev_val + next_val) / 2

            new_Section_row.append(new_val)

        new_Section_row = new_Section_row[1::]

        prev_Section_mat[-1] = prev_Section_row
        next_Section_mat[0] = next_Section_row

        new_Section_mat = np.vstack((new_Section_mat, prev_Section_mat))
        new_Section_mat = np.vstack((new_Section_mat, new_Section_row))
        new_Section_mat = np.vstack((new_Section_mat, next_Section_mat))

        new_Section_mat = new_Section_mat[1::]

        #####################################

        prev_Airfoil_mat = airfoil_Geo[:section_index + 1]
        next_Airfoil_mat = airfoil_Geo[section_index + 1::]

        new_Airfoil_mat = [['0']]

        prev_val = prev_Airfoil_mat[-1]
        next_val = next_Airfoil_mat[0]

        new_val = prev_val

        new_Airfoil_mat = np.vstack((new_Airfoil_mat, prev_Airfoil_mat))
        new_Airfoil_mat = np.vstack((new_Airfoil_mat, new_val))
        new_Airfoil_mat = np.vstack((new_Airfoil_mat, next_Airfoil_mat))

        new_Airfoil_mat = new_Airfoil_mat[1::]

        #####################################

        aero_Body.planform_Geo = new_Planform_mat.tolist()
        aero_Body.section_Geo = new_Section_mat.tolist()
        aero_Body.airfoil_Geo = new_Airfoil_mat.tolist()


        body_filename = gen_filename(user, f'{body_name}_Geo.json', save_name)

        update_aero_body(body_filename, aero_Body)

def del_body_section(aero_Body, user, save_name, section_index):

    body_name = aero_Body.name
    planform_Geo = aero_Body.planform_Geo
    section_Geo = aero_Body.section_Geo
    airfoil_Geo = aero_Body.airfoil_Geo

    index = int(section_index)

    if index == 0:
        print('Cannot delete first section')
    elif index == len(section_Geo) - 1:
        print('Cannot delete last section')
    else:

        prev_Planform_mat = planform_Geo[:index]
        next_Planform_mat = planform_Geo[index + 1::]

        del_Planform_vec = planform_Geo[index]
        del_span = del_Planform_vec[3]
        prev_Planform_mat[-1][3] = prev_Planform_mat[-1][3] + del_span

        ##################

        prev_Section_mat = section_Geo[:index]
        next_Section_mat = section_Geo[index + 1::]

        ##################

        prev_Airfoil_mat = airfoil_Geo[:index]
        next_Airfoil_mat = airfoil_Geo[index + 1::]

        ##################

        new_Planform_array = np.vstack((prev_Planform_mat, next_Planform_mat))
        new_Section_array = np.vstack((prev_Section_mat, next_Section_mat))
        new_Airfoil_array = np.vstack((prev_Airfoil_mat, next_Airfoil_mat))

        aero_Body.planform_Geo = new_Planform_array.tolist()
        aero_Body.section_Geo = new_Section_array.tolist()
        aero_Body.airfoil_Geo = new_Airfoil_array.tolist()

        body_filename = gen_filename(user, f'{body_name}_Geo.json', save_name)

        update_aero_body(body_filename, aero_Body)

def add_body(aero_Body, save_name):

    body_name = aero_Body.name
    body_orientation = aero_Body.orientation
    body_mirror = aero_Body.mirror
    body_num_hinge = aero_body.num_hinge

    body_list_csv = 'saves/{0}/Body_List.csv'.format(save_name)
    Body_headers = ['Body name', 'orientation','mirrored','num_hinge']

    Body_List = process_string_csv(body_list_csv, 0, 0)

    name_check = body_name in Body_List

    str_checc = str_check(body_name)

    if name_check == 0 and str_checc == 1:

        new_body_row = [body_name, body_orientation,body_mirror, body_num_hinge]

        new_body_list = np.vstack((Body_List, new_body_row))

        full_new_list = np.vstack((Body_headers, new_body_list))

        program_path = os.getcwd() #find path of the main file

        absolute_geo_path = program_path + '\config\Base_Save\Base_Geometry'

        base_airfoil_geo_file = r'{0}\Base_Airfoil_Geometry.csv'.format(absolute_geo_path)
        base_planform_geo_file = r'{0}\Base_Planform_Geometry.csv'.format(absolute_geo_path)
        base_section_geo_file = r'{0}\Base_Section_Geometry.csv'.format(absolute_geo_path)

        [planform_filename, section_filename, airfoil_filename] = get_filenames(save_name,body_name)

        shutil.copyfile(base_airfoil_geo_file, airfoil_filename)
        shutil.copyfile(base_planform_geo_file, planform_filename)
        shutil.copyfile(base_section_geo_file, section_filename)

        file = open(body_list_csv, 'w+', newline='')
        with file:
            write = csv.writer(file)
            write.writerows(full_new_list)

    elif name_check == 1:
        print('already body with that name')

    elif str_checc == 1:
        print('name string invalid')

    return [name_check,str_checc]

def delete_body(index, save_name):

    body_list_csv = 'saves/{0}/Body_List.csv'.format(save_name)

    Body_List = process_string_csv(body_list_csv, 0, 0)
    body_name = Body_List[index][0]

    prev_list = Body_List[:index]
    after_list = Body_List[index+1::]

    new_list = np.vstack((prev_list, after_list))

    update_file_with_mat(new_list, body_list_csv, Body_headers)

    program_path = os.getcwd()  # find path of the main file

    [Planform_filename, Section_filename, Airfoil_filename] = get_filenames(save_name, body_name)

    os.remove(Planform_filename)
    os.remove(Section_filename)
    os.remove(Airfoil_filename)

def modify_body(index, body_name, body_orientation, body_mirror, body_num_hinge, save_name):

    body_list_csv = 'saves/{0}/Body_List.csv'.format(save_name)

    Body_List = process_string_csv(body_list_csv, 0, 0)

    prev_list = Body_List[:index]
    after_list = Body_List[index + 1::]

    new_body_row = [body_name, body_orientation, body_mirror, body_num_hinge]

    new_list = np.vstack((prev_list, new_body_row, after_list))

    update_file_with_mat(new_list, body_list_csv, Body_headers)

def add_mass_item(item_name,mass,CG_FS,CG_BL,CG_WL, save_name):

    mass_csv = 'saves/{0}/Mass_Inertial_Properties.csv'.format(save_name)

    prev_inertial_list = process_string_csv(mass_csv,0,0)

    new_inertial_list = np.vstack((prev_inertial_list,[item_name,mass,CG_FS,CG_BL,CG_WL]))

    update_file_with_mat(new_inertial_list,mass_csv,inertial_headers)

def delete_mass_item(index, save_name):

    mass_csv = 'saves/{0}/Mass_Inertial_Properties.csv'.format(save_name)

    index = np.int(index)

    prev_inertial_list = process_string_csv(mass_csv, 0, 0)

    prev_list = prev_inertial_list[:index]
    after_list = prev_inertial_list[index+1::]

    new_inertial_list = np.vstack((prev_list,after_list))

    update_file_with_mat(new_inertial_list,mass_csv,inertial_headers)

def update_balance_list(index, vec, save_name):

    balance_list_csv = 'saves/{0}/Aero_Balance_List.csv'.format(save_name)

    balance_list = process_string_csv(balance_list_csv, 0, 0)


    if index < len(balance_list)-1:
        balance_list[index] = vec[0]
    elif index == len(balance_list):
        balance_list = np.vstack((balance_list, vec[0]))

    update_file_with_mat(balance_list,balance_list_csv,headers)

def make_proper_list(improper_list):

    proper_list = []
    for n in range(len(improper_list)):
        proper_list.append(improper_list[n])

    return proper_list

#need to update once file structure is fully ready
def add_save(save_name, new_aircraft_name):

    Body_List = process_string_csv(save_list_csv, 0, 0)

    name_check = save_name in Body_List

    str_checc = str_check(save_name)

    if name_check == 0 and str_checc == 1:

        new_body_row = [save_name, new_aircraft_name]

        new_body_list = np.vstack((Body_List, new_body_row))

        full_new_list = np.vstack((save_Body_headers, new_body_list))

        program_path = os.getcwd()  # find path of the main file

        file = open(save_list_csv, 'w+', newline='')
        with file:
            write = csv.writer(file)
            write.writerows(full_new_list)

        #######################################

        program_path = os.getcwd()  # find path of the main file

        print(program_path)

        new_dir_path = '{0}\saves\{1}'.format(program_path, save_name)
        os.mkdir(new_dir_path, mode=0o666)

        base_path = '{0}\config\Base_Save'.format(program_path)

        aero_bal_path = '{0}\Aero_Balance_List.csv'.format(base_path)
        body_path = '{0}\Body_List.csv'.format(base_path)
        mass_inertial_path = '{0}\Mass_Inertial_Properties.csv'.format(base_path)
        mass_load_path = '{0}\Mass_Load_Properties.csv'.format(base_path)

        geo_folder = '{0}\Base_Geometry'.format(base_path)

        new_aero_bal_path = '{0}\Aero_Balance_List.csv'.format(new_dir_path)
        new_body_path = '{0}\Body_List.csv'.format(new_dir_path)
        new_mass_inertial_path = '{0}\Mass_Inertial_Properties.csv'.format(new_dir_path)
        new_mass_load_path = '{0}\Mass_Load_Properties.csv'.format(new_dir_path)

        new_geo_folder = '{0}\Geometry_Files'.format(new_dir_path)

        shutil.copyfile(aero_bal_path, new_aero_bal_path)
        shutil.copyfile(body_path, new_body_path)
        shutil.copyfile(mass_inertial_path, new_mass_inertial_path)
        shutil.copyfile(mass_load_path, new_mass_load_path)

        shutil.copytree(geo_folder, new_geo_folder)

    elif name_check == 1:
        print('already save with that name')

    elif str_checc == 1:
        print('name string invalid')

    return [name_check, str_checc]

def delete_save(save_index):

    index = save_index

    Body_List = process_string_csv(save_list_csv, 0, 0)

    prev_list = Body_List[:index]
    after_list = Body_List[index+1::]

    new_list = np.vstack((prev_list, after_list))

    update_file_with_mat(new_list, save_list_csv, save_Body_headers)

def duplicate_save(save_index):

    num_saves = len(save_List)

    save_name = save_List[save_index][0]
    aircraft_name = save_List[save_index][1]

    new_save_name = '{0}_({1})'.format(save_name,num_saves)
    new_aircraft_name = '{0}_({1})'.format(aircraft_name,num_saves)

    new_body_row = [new_save_name, new_aircraft_name]

    new_body_list = np.vstack((save_List, new_body_row))

    full_new_list = np.vstack((save_Body_headers, new_body_list))

    program_path = os.getcwd()  # find path of the main file

    file = open(save_list_csv, 'w+', newline='')
    with file:
        write = csv.writer(file)
        write.writerows(full_new_list)

    old_save_path = '{0}\saves\{1}'.format(program_path,save_name)
    new_save_path = '{0}\saves\{1}'.format(program_path, new_save_name)

    shutil.copytree(old_save_path,new_save_path)

def modify_save(save_index, new_save_name, new_aircraft_name):

    num_saves = len(save_List)

    save_name = save_List[save_index][0]
    aircraft_name = save_List[save_index][1]

    prev_list = save_List[:save_index]
    after_list = save_List[save_index + 1::]

    new_body_row = [new_save_name, new_aircraft_name]

    new_list = np.vstack((prev_list, new_body_row, after_list))

    update_file_with_mat(new_list, save_list_csv, save_Body_headers)

    program_path = os.getcwd()  # find path of the main file
    old_save_path = '{0}\saves\{1}'.format(program_path, save_name)
    new_save_path = '{0}\saves\{1}'.format(program_path, new_save_name)

    os.rename(old_save_path,new_save_path)

def update_log(text):

    file = open('Log_File.txt', 'a')

    new_log = '{0}\n'.format(text)

    file.write(new_log)

    file.close()