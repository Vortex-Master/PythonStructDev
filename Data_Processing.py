#processing information to change its format for other functionw in the code
import csv
import numpy as np

#process a csv file down to a numpy array
def process_csv(csv_filename, start_row, start_column):

    file = open(csv_filename)
    file_opened = csv.reader(file)

    num_col = len(next(file_opened))

    processed_file = np.array(np.zeros((1, num_col)))

    for row in file_opened:
        processed_file = np.append(processed_file, [row], axis=0)

    processed_file = processed_file[1 + start_row:, start_column:]
    processed_file = processed_file.astype(np.double)

    return processed_file

#process a csv file down to a numpy list
def process_string_csv(csv_filename, start_row, start_column):

    file = open(csv_filename)
    file_opened = csv.reader(file)

    num_col = len(next(file_opened))

    processed_file = np.array(np.zeros((1, num_col)))

    for row in file_opened:
        processed_file = np.append(processed_file, [row], axis=0)

    processed_file = processed_file[1 + start_row:, start_column:]

    return processed_file

def process_airfoil_pts(airfoil_pts):

    prev_x = 0
    row_x = .5 #initialize it to something we know wont work for now

    n = 0
    while prev_x != row_x:

        prev_x = row_x

        row_x = airfoil_pts[n][0]

        n += 1

    new_airfoil_pts = airfoil_pts[0:n-1, :]

    return new_airfoil_pts

#generate an array for new airfoil to be added to database
def process_pasted_airfoil_data(airfoil_name, str_airfoil_pts):

    start_array = [[airfoil_name, ''],
                   ['ac', ''],
                   ['x', 'z']]

    full_array = [['','']]

    real_vals = ['0']

    working_val = ''

    for i in range(len(str_airfoil_pts)):

        val = str_airfoil_pts[i]

        if val == ' ':

            if working_val != '':

                new_val = np.double(working_val)

                real_vals.append(new_val)

                working_val = ''

        else:

            working_val += val

    # take off the initail zero and add zero to the end
    del real_vals[0]
    real_vals.append(0)

    x_vals = real_vals[::2]
    z_vals = real_vals[::-2]
    z_vals = z_vals[::-1]

    x_vals[0] = 1
    x_vals[-1] = 1
    z_vals[0] = 0
    z_vals[-1] = 0

    pts_mat = np.vstack((x_vals,z_vals)).T

    zerozero_diff = abs(pts_mat[:,0])
    zero_pt_min_val = min(zerozero_diff)
    zero_pt_min_ind = np.where(zerozero_diff == zero_pt_min_val)[0][0]
    pts_mat[zero_pt_min_ind] = [0,0]

    print(pts_mat)

    full_array = np.vstack((start_array,pts_mat))

    return [full_array, pts_mat]

