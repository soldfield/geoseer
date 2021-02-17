"""
Scripts to import data from an ascii formatted file containing multiple polyline geometries.

The file should have no header, with data starting on the first line.

Each polyline should be listed as a tab-deliminated string of coordinates, with either two or three coordinates per line.
These columns should be X and Y (and optionally Z).

Despite handling 3D coordinates, the code assumes that the dataset may be safely projected to 2D.
This is will work best with datasets where the surface upon which polylines/lineations are interpreted is approaching
a planar surface. For example a cliff face that is more extensive laterally and vertically but has limited offset
 perpendicular to the cliff-face (along the line of sight of the observer).

"""
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import MultiLineString
# from shapely.geometry import LineString, asLineString, asMultiLineString
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
    Makes a shapely MultiLineString object from list of coordinates
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
    "Plots all lines extracted from file on X and Y axes."
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

input_file = "C:\\Users\\simold\\Documents\\git\\DFNcompare\\data\\tala_cc_fracs\\measurements.poly"
file_path = os.path.normpath(input_file)
coord_list = extract_coordinates_from_xyz_file(file_path)
all_lines = make_linestring_list(coord_list)
bounding_rectangle = calculate_bounding_rectangle(all_lines)
# plot_all_fracs_xy_plane(all_lines)
angle = bounding_box_long_axis(bounding_rectangle)

#%%

def plot_bounding_box(reproj_bound_box):
    x, y = reproj_bound_box.exterior.xy
    plt.plot(x,y)
    plt.show()

def reproject_to_local_crs(all_lines, angle, origin):
    # rotate original interpreted lines
    all_lines_reproj = rotate(all_lines, angle, origin=(560420,6323600), use_radians=False)

    for ln in all_lines_reproj:
        x = np.array([c[0] for c in ln.coords])
        x = x - 560400
        y = np.array([c[1] for c in ln.coords])
        y = y - 6323600
        z = np.array([c[2] for c in ln.coords])
        plt.plot(x,y)

    plt.ylabel("Lateral distance from lakeside/West (m)")
    plt.xlabel("Lateral distance from left/North (m)")
    plt.show()
    return all_lines_reproj

# TODO: automate identification of the origin
origin=(560420,6323600)
reproj_bound_box = rotate(bounding_rectangle, angle, origin=(560420,6323600), use_radians=False)
all_lines_reproj = reproject_to_local_crs(all_lines, angle, origin)

#%%
for ln in all_lines_reproj:
    x = np.array([c[0] for c in ln.coords])
    x = x - 560400
    y = np.array([c[1] for c in ln.coords])
    y = y - 6323600
    z = np.array([c[2] for c in ln.coords])

    plt.plot(x,z)

plt.ylabel("Vertical distance from base (m)")
plt.xlabel("Lateral distance from left/North (m)")
plt.show()

# left should be north, this has been confirmed with physical data

#%%

# plot wall section
for ln in all_lines_reproj:
    x = np.array([c[0] for c in ln.coords])
    x = x - 560400
    y = np.array([c[1] for c in ln.coords])
    y = y - 6323600
    z = np.array([c[2] for c in ln.coords])

    plt.plot(x,z)

plt.xlim(0, 50)
plt.ylabel("Vertical distance from base (m)")
plt.xlabel("Lateral distance from left/North (m)")
plt.show()

# project for a section of wall and export that data for use in fraqpak

#%%

def search_fractures_x_window(window_x_min, window_x_max, all_lines_reproj):
    """
    Script to identify fractures that exist between two defined x values.
    All other coordinates ignored.
    Returns a list of shapely linestrings
    """

    active_lines = []
    for line in all_lines_reproj:
        x_list, y_list = line.xy
        x_list = [i -560400 for i in x_list]
        line_x_max = np.max(x_list)
        line_x_min = np.min(x_list)
        # if line exists in search window do something
        if line_x_min <= window_x_max and line_x_max >= window_x_min:
            active_lines.append(line)
    return active_lines


# print(active_lines)
active_line = search_fractures_x_window(window_x_min, window_x_max, all_lines_reproj)

#%%

# Plot active lines

for ln in active_line:
    x = np.array([c[0] for c in ln.coords])
    x = x - 560400
    y = np.array([c[1] for c in ln.coords])
    y = y - 6323600
    z = np.array([c[2] for c in ln.coords])

    plt.plot(x,z)

#plt.xlim(0, 50)
plt.ylabel("Vertical distance from base (m)")
plt.xlabel("Lateral distance from left/North (m)")
plt.show()


#%%
# pipe to fracpaqpy

