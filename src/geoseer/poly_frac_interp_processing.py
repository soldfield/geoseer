"""
Scripts to import data from an ascii formatted file containing multiple polyline geometries.

The file should have no header, with data starting on the first line.

Each polyline should be listed as a tab-deliminated string of coordinates,
with either two or three coordinates per line.
These columns should be X and Y (and optionally Z).

Despite handling 3D coordinates, the code assumes that the dataset may be safely projected to 2D.
This is will work best with datasets where the surface upon which polylines/lineations are interpreted is approaching
a planar surface. For example a cliff face that is more extensive laterally and vertically but has limited offset
 perpendicular to the cliff-face (along the line of sight of the observer).

"""
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import MultiLineString
from shapely.affinity import rotate
import os

def extract_coordinates_from_xyz_file(file_path):
    with open(file_path) as f_in:
        coord_list = f_in.readlines()

    # strip newlines
    coord_list = [i.strip() for i in coord_list]
    coord_list = [i.split() for i in coord_list]
    coord_list = [np.array(i).astype(float) for i in coord_list]
    return coord_list

def make_linestring_list(coord_list):
    """
    Makes a shapely MultiLineString object from list of coordinates.
    :param coord_list:
    :return:
    """
    # converts list of coordinates into a list of lines, each stored as a numpy array
    line_list = []
    line = []
    for i in range(len(coord_list)):
        if len(coord_list[i]) != 0:
            line.append(coord_list[i])
        else:
            line = np.stack(line)
            line_list.append(line)
            line = []

    # Make collection of lines
    all_lines = MultiLineString(line_list)
    return all_lines

def calculate_bounding_rectangle(all_lines):
    # plot bounding box to fractures
    bounding_rectangle = all_lines.minimum_rotated_rectangle
    x, y = bounding_rectangle.exterior.xy
    plt.plot(x,y)
    plt.show()
    return bounding_rectangle

def plot_all_fracs_xy_plane(all_lines):
    """
    Plots all lines extracted from file on X and Y axes.
    """
    #plot individual fractures
    for i in all_lines:
        x, y = i.xy
        plt.plot(x,y)
    plt.show()

def bounding_box_long_axis(bounding_rectangle):
    # find trend
    x, y = bounding_rectangle.exterior.xy

    x1, x2, x3, x4 = x[0:4]
    y1, y2, y3, y4 = y[0:4]

    xdiff = np.abs(np.mean([x1, x4])-np.mean([x2, x3]))
    ydiff = np.abs(np.mean([y1, y4])-np.mean([y2, y3]))

    angle = np.rad2deg(np.arctan(ydiff/xdiff))
    print(angle)
    # positive angles will become counter-clockwise rotation
    return angle

def plot_bounding_box(reproj_bound_box):
    x, y = reproj_bound_box.exterior.xy
    plt.plot(x,y)
    plt.show()

def reproject_to_local_crs(all_lines, angle, origin):
    """
    Reprojects coordinates to the provided origin and rotation angle.
    Returns the reprojected coordinates in a shapely MultiLineString.
    """
    # rotate original interpreted lines
    all_lines_reproj = rotate(all_lines, angle, origin, use_radians=False)

    output_lns = []

    for i in range(len(all_lines_reproj)):
        ln = all_lines_reproj[i]
        x = np.array([c[0] for c in ln.coords])
        x = x - origin[0]
        y = np.array([c[1] for c in ln.coords])
        y = y - origin[1]
        z = np.array([c[2] for c in ln.coords])
        if len(origin) == 3:
            z = z - origin[2]
        plt.plot(x, y)

        #
        reproj_ln = list(zip(x, y, z))
        output_lns.append(reproj_ln)


    output_lns = MultiLineString(output_lns)

    plt.ylabel("Lateral distance from lakeside/West (m)")
    plt.xlabel("Lateral distance from left/North (m)")
    plt.show()
    return output_lns

def plot_xz_linelist(line_list):
    """
    Plots x and z of all lines in list.
    Assumes three coords with x and z in first and third column, respectively.
    :param line_list:
    """
    for ln in line_list:
        x = np.array([c[0] for c in ln.coords])
        #y = np.array([c[1] for c in ln.coords])
        z = np.array([c[2] for c in ln.coords])

        plt.plot(x, z)

    # plt.xlim(0, 50) # limits section of wall presented
    plt.ylabel("Vertical distance from base (m)")
    plt.xlabel("Lateral distance from left/North (m)")
    plt.show()


def search_fractures_x_window(window_x_min, window_x_max, all_lines_reproj):
    """
    Script to identify fractures that exist between two defined x values.
    All other coordinates ignored.

    :param window_x_min: start of x window
    :param window_x_max: end of x window
    :return active_lines: A list of shapely linestrings
    """

    active_lines = []
    for line in all_lines_reproj:
        x_list, y_list = line.xy
        line_x_max = np.max(x_list)
        line_x_min = np.min(x_list)
        # if line exists in search window do something
        if line_x_min <= window_x_max and line_x_max >= window_x_min:
            active_lines.append(line)

    active_lines = MultiLineString(active_lines)
    return active_lines

def full_frac_interp_process_example(input_file, origin, plot=False):
    """
    Function that calls each of the functions above to convert a poly file
    into a list of polylines each of which represent a fracture and output a 
    figure in the XZ plane.
    
    :param input_polyfile: file from cloud compare in our use case
    :return all_lines_reproj: list of polylines
    """
    
    
    file_path = os.path.normpath(input_file)
    coord_list = extract_coordinates_from_xyz_file(file_path)
    all_lines = make_linestring_list(coord_list)
    bounding_rectangle = calculate_bounding_rectangle(all_lines)

    # plot_all_fracs_xy_plane(all_lines)
    angle = bounding_box_long_axis(bounding_rectangle)

    # TODO: automate origin identfication
    # origin = (560420,6323600,0) # for Rørdal, note X, Y, Z, Z can be changed to match other figures
    
    # reproj_bound_box = rotate(bounding_rectangle, angle, origin, use_radians=False)
    all_lines_reproj = reproject_to_local_crs(all_lines, angle, origin)  # includes plot

    if plot == True:
        plot_xz_linelist(all_lines_reproj, origin)

    return all_lines_reproj


def xy_to_fracpaq(line_list: list, out_path: str) -> None:
    """
    Converts a shapely multilinestring object to a FracPaQ file, extracting data to a two-dimensional plane
    oriented to intersect the x and y axes (xy-plane).

    :param line_list: shapely Multilinestring object of lines defined by two or three coordinates, XY(Z)
    :param out_path: string of file path, name and any extensions (e.g. '.txt')
    :return None: generates text file in fracpaq format at out_path
    """
    string_out = ""

    for i in range(len(line_list)):
        # retrieve coordinates for line
        ln_coords = np.array(line_list[i].coords)
        i_arr = ln_coords[:, 0]
        j_arr = ln_coords[:, 1]
        plt.plot(i_arr, j_arr)

        for i in range(len(i_arr)):
            string_out += str(i_arr[i]) + "\t" + str(j_arr[i]) + "\t"
            # print(i, len(i_arr))
            if i == len(i_arr)-1:
                # print(i)
                string_out += "\n"
                # list_string_lines_out.append(string_line_out)
    # write to output file
    with open(out_path, "w") as f_out:
        f_out.write(string_out)

    plt.ylabel("Northing (m)")
    plt.xlabel("Easting (m)")
    plt.show()

    plot_file_name = out_path[:-4]+".png"
    plt.savefig(plot_file_name)

def xz_to_fracpaq(line_list: list, out_fig_name: str) -> None:
    """
    Converts a shapely multilinestring object to a FracPaQ file, extracting data to a two-dimensional plane
    oriented to intersect the x and y axes (xy-plane).

    :param line_list: shapely Multilinestring object of lines defined by two or three coordinates, XY(Z)
    :param out_path: string of file path, name and any extensions (e.g. '.txt')
    :return None: generates text file in fracpaq format at out_path
    """
    string_out = ""

    for i in range(len(line_list)):
        # retrieve coordinates for line
        ln_coords = np.array(line_list[i].coords)
        i_arr = ln_coords[:, 0]
        j_arr = ln_coords[:, 2]
        plt.plot(i_arr, j_arr)

        for i in range(len(i_arr)):
            string_out += str(i_arr[i]) + "\t" + str(j_arr[i]) + "\t"
            # print(i, len(i_arr))
            if i == len(i_arr)-1:
                # print(i)
                string_out += "\n"
                # list_string_lines_out.append(string_line_out)
    # write to output file
    with open(out_path, "w") as f_out:
        f_out.write(string_out)

    plt.ylabel("Vertical distance (m)")
    plt.xlabel("Horizontal distance (m)")

    plt.savefig(out_fig_name)
    plt.show()

def find_min_x(line_list):
    min_x_ls = []

    for line in line_list:
        min_x_ls.append(line.bounds[0])

    min_x = np.min(min_x_ls)
    return min_x

# input_file = "C:\\Users\\simold\\Documents\\git\\DFNcompare\\data\\tala_cc_fracs\\measurements.poly"
input_file = "C:\\Users\\shl459\\Desktop\\DTU_20201021\\measurements.poly"

origin = [560315,6323600,0]
all_lines_reprojected = full_frac_interp_process_example(input_file, origin, False)

#%%

# Run lines 1-275 once
# Define window size and wall length 
# Then run below once
# Will output to working directory (where this code is saved)

# Looped output of sample windows

window_size = int(10)
quarry_wall_length = int(1000)
number_of_windows = 2 #int(quarry_wall_length / window_size)
out_directory = os.getcwd()

line_list = all_lines_reprojected

for i in range(0,number_of_windows):
    plt.clf()

    min_x = find_min_x(line_list)

    window_x_min = (window_size * i) + min_x
    window_x_max = (window_size * (i +1)) + min_x

    name = "fpq_input_window_" + str(window_x_min) + "_to_" + str(window_x_max)
    out_file_name = name +'.txt'
    out_path = os.path.join(out_directory, out_file_name)

    fx_lines = search_fractures_x_window(window_x_min, window_x_max, line_list)

    out_fig_name = name +'.png'
    xz_to_fracpaq(fx_lines, out_fig_name)



#%%