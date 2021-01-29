

"""
   Copyright 2020 Simon Oldfield

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""


from shapely.geometry import LineString
import numpy as np
from matplotlib import pyplot as plt
import os


#%%

def seg_len(pt1, pt2):
    x1, y1 = pt1[0], pt1[1]
    x2, y2 = pt2[0], pt2[1]

    seg_len = np.sqrt(np.power(np.abs(x1 - x2), 2) + np.power(np.abs(y1 - y2), 2))

    return seg_len

def running_len(coord_list):
    i = 0
    len = 0
    while i < len(coord_list)-1:
        x1, y1 = coord_list[i]
        x2, y2 = coord_list[i+1]
        i += 1
        len += seg_len((x1, y1), (x2, y2))
        print(i)
        yield len


#%%

#plt.clf()

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.scatter(y*1000,x*1000)
ax.ticklabel_format(useOffset=False)

plt.show()


#%%
plt.clf()
a = LineString([(0,0),(1,5),(5,10)])

i, j = a.xy

plt.plot(a.coords)
pt = a.interpolate(1).xy
plt.scatter(pt[0], pt[1])
plt.show()


#%%


class SurveyMarker:
    """holder for location data of survey markers"""
    def __init__(self, id, x,y):
        self.id = id
        self.x = x
        self.y = y

class SurveyLine:
    def __init__(self, line_name):
        self.line_name = line_name

    def add_trace_id(self, start, end):
        self.start_trace = start
        self.end_trace = end

    def add_markers_byfile(self, file_path):

        self.marker_list = []
        ln_name = os.path.basename(file_path).split('_')[0]

        with open(src, "r") as data:
            db = data.read()

        db = db.split('\n')
        for i in range(len(db)):
            db[i] = db[i].split()



        np.array(db).astype('float')

#%%



src_path = "O:\Public\DFN\Validation\Rordal\gpr_survey"
file_ls = os.listdir(src_path)[-7:-2]

#%%

marker_dict = {}
for file in file_ls:

    src = os.path.join(src_path, file)

    line_name = os.path.basename(src).split('_')[0]
    line_ltr = line_name[-1]

    with open(src, "r") as data:
        db = data.read()

    db = db.split('\n')
    for i in range(len(db)):
        db[i] = db[i].split()
        db[i] = db[i]
        marker_dict[line_ltr+str(i+1)] = db[i]


print(marker_dict)


#%%

trc_loc_file = "MarkerTraceLocations.txt"
src = os.path.join(src_path, trc_loc_file)

with open(src, "r") as data:
    db = data.read()

db = db.split('\n')
for i in range(len(db)):
    db[i] = db[i].split('\t')
    db[i][2] = float(db[i][2])


#%%

# append xy coord to each marker trace location

mkr_trc_locs = db

line_list = np.unique([i[0] for i in mkr_trc_locs])

# append coordinates
for i in range(len(mkr_trc_locs)):
    coords = marker_dict[mkr_trc_locs[i][1]]
    mkr_trc_locs[i].append(coords)

#%%

ln_name = 'Line_15'

def line_point_extraction(ln_name, mkr_trc_locs):
    """

    :param ln_name:
    :param mkr_trc_locs:
    :return:
    """
    ln_pt_ls = []
    for i in range(len(mkr_trc_locs)):
        if mkr_trc_locs[i][0]==ln_name:
            ln_pt_ls.append(i)

    line_pts = None
    line_pts = [mkr_trc_locs[i] for i in ln_pt_ls]
    # Assumes that data is in the desired order

    #option to sort if needed:
    #sorted_line_pts = sorted(line_pts, key=lambda x: x[2])

    return line_pts

line_pts = line_point_extraction(ln_name, mkr_trc_locs)

#%%

def line_seg_point_interpolation(line_pts):
    """
    Fucntion to interpolate trace coordinates from a dictionary containing
    line name, marker names and marker coordinates.

    :param line_pts:
    :return:
    """
    line_seg_pts = []

    for i in range(len(line_pts)):
        # add known point

        line_seg_pts.append([float(num) for num in line_pts[i][3]])

        if i == int(len(line_pts)-1):
            continue
        else:
            # extract point coordinates
            x1, y1 = [float(num) for num in line_pts[i][3]]
            x2, y2 = [float(num) for num in line_pts[i+1][3]]

            # calc segment lateral diffs
            x_diff = x2 - x1
            y_diff = y2 - y1

            trace_num = int(line_pts[i + 1][2] - line_pts[i][2])

            x_spac = x_diff / trace_num
            y_spac = y_diff / trace_num

            # loop for addition of points
            for j in range(trace_num-1):
                x_n = x1 + (x_spac * j)
                y_n = y1 + (y_spac * j)
                line_seg_pts.append([x_n, y_n])

    print(line_seg_pts)
    return line_seg_pts

line_seg_pts = line_seg_point_interpolation(line_pts)


# plot result
# check lines that were shot in reverse (i.e. A5 to A1 etc.)

# export text file

#%%

x_series = [pt[0] for pt in line_seg_pts]
y_series = [pt[1] for pt in line_seg_pts]


fig, ax = plt.subplots()

ax.scatter(x_series, y_series)
#ax.set_aspect('equal')

#plt.plot([('9.997917', '57.050412'),('9.997823', '57.050406'),('9.997721', '57.050406')])

ax.grid(True)
plt.xlim(9.99771, 9.99792)
fig.tight_layout()

plt.show()

#%%

line_pts_dict = {}

for ln_name in line_list:
    line_pts = line_point_extraction(ln_name, mkr_trc_locs)
    line_seg_pts = line_seg_point_interpolation(line_pts)
    line_pts_dict[ln_name] = line_seg_pts

#%%


for key in line_pts_dict.keys():

    #key = 'Line_13'
    data = line_pts_dict[key]


    src_path = "O:\Public\DFN\Validation\Rordal\gpr_survey"
    fn = key+"_navdat.txt"
    f_path = os.path.join(src_path, fn)

    np.savetxt(f_path, np.array(data))
