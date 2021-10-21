
import numpy as np
from shapely.geometry import MultiLineString
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


input_file = "C:\\Users\\simold\\Documents\\git\\DFNcompare\\data\\tala_cc_fracs\\measurements.poly"


coord_list = extract_coordinates_from_xyz_file(input_file)
lines = make_linestring_list(coord_list)


#%%

f = open("test.txt", "w")


frac_count = len(lines)

# write format header block

f.writelines("BEGIN Format\n")
f.writelines("Format = Ascii\n")
f.writelines("Scale = 100.\n")
f.writelines(f"NumFractures = {frac_count}\n")
f.writelines("END Format\n")

# write properties block
f.writelines("BEGIN Properties\n")
f.writelines("Prop1 = (Real * 4)      \"Transmissivity\"\n")
f.writelines("Prop2 = (Real * 4)      \"Storativity   \"\n")
f.writelines("Prop3 = (Real * 4)      \"Frac Thickness\"\n")
f.writelines("END Properties\n")

# write fracture data block
f.writelines("BEGIN Fracture\n")

frac_records = []
frac_num_i = 0

for li in range(len(lines[0:10])):
    frac_coords = list(lines[l].coords)
    frac_num_li = li+1
    vertex_count = len(frac_coords)
    frac_set = 0
    f.writelines(f"{frac_num_li:6}{vertex_count:6}{frac_set:6}{1:6}{1:6}{1:6}\n")
    for vi in range(len(frac_coords)):
        f.writelines(f"{vi+1:6}{frac_coords[vi][0]:16.6}{frac_coords[vi][1]:16.6}{frac_coords[vi][2]:16.6}\n")

f.writelines("END Fracture\n")

f.close()

#%%


