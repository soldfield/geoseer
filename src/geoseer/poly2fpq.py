"""
A script to convert poly files of polylines into a readable formwat for fracpaq.
"""

import os
folder_path = os.getcwd()
input_file = "measurements.poly"

def poly2fpq(input_file):
    """
    Function accepts .poly files and converts them to a format acceptable to FracPaQ (Healy et al. 2017).
    The first two columns of coordinates will be treated as X and Y, respectively. Any third field for Z is discarded.
    :param input_file: file name of the .poly file
    :return: Generates a text file in the current working directory, named as the input file with suffix "_poly2fpq"
    """
    with open(input_file) as f_in:
        coord_list = f_in.readlines()
        #print(coord_list)

    # strip newlines
    coord_list = [i.strip() for i in coord_list]
    coord_list = [i.split() for i in coord_list]
    #test = [[float(j) for j in test[i]] for i in test]

    # format each record into a list of points: x1 y1 ... xn yn
    for i in range(len(coord_list)):
        if len(coord_list[i]) > 0:
            coord_list[i] = coord_list[i][0] + " " + coord_list[i][1] + " "
        else:
            coord_list[i] = "\n"

    # compile each list of points into a list of lines
    db = ""
    for i in range(len(coord_list)):
        db = db + coord_list[i]

    # write to output file
    output_file_name = input_file.strip('.poly')+"_poly2fpq.txt"
    with open(output_file_name, "w") as f_out:
        f_out.write(db)

def xz_to_fracpaq(line_list):
    """
    Converts a shapely multilinestring object to a FracPaQ file, extracting data to a two-dimensional plane
    oriented to intersect the x and z axes (xz-plane).

    :param line_list: shapely Multilinestring object of lines defined by three coordinates, XYZ
    :param out_path: string of file path, name and any extensions (e.g. '.txt')
    :return None: generates text file in frcpaq format at out_path
    """

    try:
        line_list.has_z
    except ValueError:
        print("Coordinates missing third dimension, no z value")

    for i in range(len(line_list)):
        # retrieve coordinates for line
        ln_coords = np.array(line_list[i].coords)
        x_arr = ln_coords[:,0]
        z_arr = ln_coords[:,2]
        plt.plot(x_arr, z_arr)

    plt.ylabel("Vertical distance from base (m)")
    plt.xlabel("Lateral distance from left/North (m)")
    plt.show()

line_list = active_lines
xz_to_fracpaq(line_list)

#%%

def xy_to_fracpaq(line_list: list, out_path: str) -> None:
    """
    Converts a shapely multilinestring object to a FracPaQ file, extracting data to a two-dimensional plane
    oriented to intersect the x and y axes (xy-plane).

    :param line_list: shapely Multilinestring object of lines defined by two or three coordinates, XY(Z)
    :param out_path: string of file path, name and any extensions (e.g. '.txt')
    :return None: generates text file in frcpaq format at out_path
    """
    string_out = ""

    for i in range(len(line_list)):
        # retrieve coordinates for line
        ln_coords = np.array(line_list[i].coords)
        i_arr = ln_coords[:, 0]
        j_arr = ln_coords[:, 1]
        plt.plot(i_arr, j_arr)

        for i in range(len(i_arr)):
            string_out += str(i_arr[i]) + " " + str(j_arr[i]) + " "
            #print(i, len(i_arr))
            if i == len(i_arr)-1:
                #print(i)
                string_out += "\n"
                # list_string_lines_out.append(string_line_out)
    # write to output file
    with open(out_path, "w") as f_out:
        f_out.write(string_out)

    plt.ylabel("Northing (m)")
    plt.xlabel("Easting (m)")
    plt.show()

out_dir = os.getcwd()
out_file_name = 'test.txt'
out_path = os.path.join(out_dir, out_file_name)

line_list = active_lines
xy_to_fracpaq(line_list, out_path)

#%%

def coord_arrays_to_fpq(i_arr, j_arr, out_path):
    """
    Outputs a text file in fracpaq format from two lists of coordinates.

    :param i_arr: array of coordinates for horizontal axis i.e. x
    :param j_arr: array of coordinates for vertical axis i.e. y or z
    """
    # format each record into a list of points: x1 y1 ... xn yn
    coord_list = []

    for i in range(len(i_arr)):
        if len(i_arr[i]) > 0:
            coord_list.append(i_arr[i] + " " + j_arr[i] + " ")
        else:
            coord_list.append("\n")

    # compile each list of points into a list of lines
    db = ""
    for i in range(len(coord_list)):
        db = db + coord_list[i]

    # write to output file
    with open(out_path, "w") as f_out:
        f_out.write(db)

out_dir = os.getcwd()
out_file_name = 'test.txt'
out_path = os.path.join(out_dir, out_file_name)

i_arr = x_arr
j_arr = z_arr

coord_arrays_to_fpq(i_arr, j_arr, out_path)


#%%


input_file = "/data/tala_cc_fracs/measurements.poly"


poly2fpq(input_file)

