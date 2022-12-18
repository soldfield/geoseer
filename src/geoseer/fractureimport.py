"""Fracture import module"""

import os
from tkinter import filedialog as fd

import numpy as np
import pandas as pd
import shapely.geometry as sg
import shapely.affinity as sa
import shapely.ops as so
import matplotlib.pyplot as plt

from poly_frac_interp_processing import *


class ImportFracFile:
    """Class for initial import of fracture files."""

    def __init__(self, file_path=None):
        if file_path == None:
            self.file_path = os.path.normpath(fd.askopenfilename())
        else:
            self.file_path = file_path

        self.file_name = os.path.basename(file_path)
        assert os.path.isfile(self.file_path), "Source is not recognised as a file"
        assert (
            self.file_name[-5:] == ".poly"
        ), "Source file name does not end in '.poly'"

        self.import_data = self.read_coordinates_from_xyz_file()
        self.make_linestring_list()
        self.get_data_long_axis()

    def read_coordinates_from_xyz_file(self):
        with open(self.file_path) as f_in:
            coord_list = f_in.readlines()
            # strip newlines
        coord_list = [i.strip() for i in coord_list]
        coord_list = [i.split() for i in coord_list]
        coord_list = [np.array(i).astype(float) for i in coord_list]
        self.coord_list = coord_list

    def make_linestring_list(self):
        """
        Makes a shapely MultiLineString object from list of coordinates.
        :param coord_list:
        :return:
        """
        # converts list of coordinates into a list of lines, each stored as a numpy array
        line_list = []
        line = []
        for i in range(len(self.coord_list)):
            if len(self.coord_list[i]) != 0:
                line.append(self.coord_list[i])
            else:
                line = np.stack(line)
                line_list.append(line)
                line = []

        # Make collection of lines
        self.original_lines = sg.MultiLineString(line_list)

    def get_data_long_axis(self):
        min_bounding_rectangle = self.original_lines.minimum_rotated_rectangle
        x, y = min_bounding_rectangle.exterior.xy
        x1, x2, x3, x4 = x[0:4]
        y1, y2, y3, y4 = y[0:4]

        xdiff = np.abs(np.mean([x1, x4]) - np.mean([x2, x3]))
        ydiff = np.abs(np.mean([y1, y4]) - np.mean([y2, y3]))

        self.long_axis_angle = np.rad2deg(np.arctan(ydiff / xdiff))

    def reproject_to_local_grid(self):
        """
        Reprojects coordinates to the provided origin and rotation angle.
        Returns the reprojected coordinates in a shapely MultiLineString.
        """
        angle = self.long_axis_angle
        original_lines = self.original_lines
        self.original_extent = original_lines.bounds

        # rotate original interpreted lines
        rot_lns = sa.rotate(original_lines, angle, origin="center", use_radians=False)

        minx, miny, maxx, maxy = rot_lns.bounds
        self.rotated_extent = [minx, miny, maxx, maxy]

        # translate: geometry, x offset, y offset, z offset
        rot_lns = sa.translate(rot_lns, -minx, -miny, 0)

        self.translated_extent = rot_lns.bounds

        self.reprojected_lines = rot_lns

        # plt.ylabel("Lateral distance from lakeside/West (m)")
        # plt.xlabel("Lateral distance from left/North (m)")
        # plt.show()
        return rot_lns

    def plot_xz_linelist(self):
        """
        Plots x and z of all lines in list.
        Assumes three coords with x and z in first and third column, respectively.
        :param line_list:
        """

        line_list = self.reprojected_lines

        for ln in line_list:
            x = np.array([c[0] for c in ln.coords])
            # y = np.array([c[1] for c in ln.coords])
            z = np.array([c[2] for c in ln.coords])

            plt.plot(x, z)

        # plt.xlim(0, 50) # limits section of wall presented
        plt.ylabel("Vertical distance from base (m)")
        plt.xlabel("Lateral distance from left/North (m)")
        plt.show()

    def search_fractures_x_window(window_x_min, window_x_max):
        """
        Script to identify fractures that exist between two defined x values.
        All other coordinates ignored.

        :param window_x_min: start of x window
        :param window_x_max: end of x window
        :return active_lines: A list of shapely linestrings
        """

        all_lines_reproj = self.reprojected_lines

        active_lines = []
        for line in all_lines_reproj:
            x_list, y_list = line.xy
            line_x_max = np.max(x_list)
            line_x_min = np.min(x_list)
            # if line exists in search window do something
            if line_x_min <= window_x_max and line_x_max >= window_x_min:
                active_lines.append(line)

        active_lines = sg.MultiLineString(active_lines)
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

        # TODO: automate origin identification
        # origin = (560420,6323600,0) # for Rørdal, note X, Y, Z, Z can be changed to match other figures

        # reproj_bound_box = rotate(bounding_rectangle, angle, origin, use_radians=False)
        all_lines_reproj = reproject_to_local_crs(
            all_lines, angle, origin
        )  # includes plot

        if plot == True:
            plot_xz_linelist(all_lines_reproj, origin)

        self.all_lines_reproj = all_lines_reproj

        return all_lines_reproj

    def dip_from_ij_arr(i_arr, j_arr) -> float:
        """
        Returns a dip angle in dagrees from the start and end point of a line defined by a list of i and j coordinates
        :param i_arr: list of i coordinates
        :param j_arr: list of j coordinates
        :return: float of dip in degrees
        """

        p1_x, p2_x = i_arr[0], i_arr[-1]
        p1_y, p2_y = j_arr[0], j_arr[-1]

        max_x = max(p1_x, p2_x)
        min_x = min(p1_x, p2_x)

        max_y = max(p1_y, p2_y)
        min_y = min(p1_y, p2_y)

        m = (max_y - min_y) / (max_x - min_x)

        dip = math.degrees(math.atan(m))

        return dip

    def rgba_col_val(value, cmap_name="Spectral", scale_min=0.0, scale_max=90.0):
        cmap = mpl.cm.get_cmap(cmap_name)  # insert colormap name here to change
        norm = mpl.colors.Normalize(vmin=scale_min, vmax=scale_max)
        rgb_val = cmap(norm(value))
        return rgb_val

    def print_colorbar(cmap_name="viridis", scale_min=0.0, scale_max=90.0) -> None:
        cmap = mpl.cm.get_cmap(cmap_name)  # insert colormap name here to change
        norm = mpl.colors.Normalize(vmin=scale_min, vmax=scale_max)
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation="horizontal",
            label="Dip (degrees)",
        )
        plt.show()

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
                if i == len(i_arr) - 1:
                    # print(i)
                    string_out += "\n"
                    # list_string_lines_out.append(string_line_out)
        # write to output file
        with open(out_path, "w") as f_out:
            f_out.write(string_out)

        plt.ylabel("Northing (m)")
        plt.xlabel("Easting (m)")
        plt.show()

        plot_file_name = out_path[:-4] + ".png"
        plt.savefig(plot_file_name)

    def xz_to_fracpaq(line_list: list, out_fig_name: str, xlim: list) -> None:
        """
        Converts a shapely multilinestring object to a FracPaQ file, extracting data to a two-dimensional plane
        oriented to intersect the x and y axes (xy-plane).

        :param line_list: shapely Multilinestring object of lines defined by two or three coordinates, XY(Z)
        :param out_path: string of file path, name and any extensions (e.g. '.txt')
        :return None: generates text file in fracpaq format at out_path
        """
        string_out = ""

        fig, ax = plt.subplots()
        # fig.figure(figsize=(8,3), dpi=300) # change figure size and resolution

        for i in range(len(line_list)):
            # retrieve coordinates for line
            ln_coords = np.array(line_list[i].coords)
            i_arr = ln_coords[:, 0]
            j_arr = ln_coords[:, 2]

            dip_ij = dip_from_ij_arr(i_arr, j_arr)
            dip_rgba = rgba_col_val(
                dip_ij, cmap_name="viridis", scale_min=0.0, scale_max=90.0
            )

            ax.plot(i_arr, j_arr, color=dip_rgba)

            for i in range(len(i_arr)):
                string_out += str(i_arr[i]) + "\t" + str(j_arr[i]) + "\t"
                # print(i, len(i_arr))
                if i == len(i_arr) - 1:
                    # print(i)
                    string_out += "\n"
                    # list_string_lines_out.append(string_line_out)
        # write to output file
        with open(out_path, "w") as f_out:
            f_out.write(string_out)

        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylabel("Vertical distance (m)")
        ax.set_xlabel("Horizontal distance (m)")
        ax.set_aspect("equal")

        fig.savefig(out_fig_name)
        plt.show()

    def find_min_x(line_list, round_down=True):
        min_x_ls = []

        for line in line_list:
            min_x_ls.append(line.bounds[0])

        min_x = np.min(min_x_ls)

        if round_down == True:
            min_x = math.floor(min_x)

        return min_x


# # obj = ImportFracFile()
# # print(obj.long_axis_angle)
