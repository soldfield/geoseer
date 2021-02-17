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


#%%


input_file = "/data/tala_cc_fracs/measurements.poly"


poly2fpq(input_file)

