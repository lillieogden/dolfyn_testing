# To get started first import the DOLfYN ADV advanced programming
# interface (API):
import dolfyn.adv.api as avm

# Import matplotlib tools for plotting the data:
from matplotlib import pyplot as plt
import matplotlib.dates as dt
import numpy as np

# The file to load:
fname = '/Users/lillie/turbulence_data/raw_data/TTM_NREL03_May2015.VEC'

# This is the vector from the ADV head to the body frame, in meters,
# in the ADV coordinate system.
body2head_vec = np.array([9.75, 2, -5.75]) * 0.0254


# This is the orientation matrix of the ADV head relative to the body.
# In this case the head was aligned with the body, so it is the
# identity matrix:
body2head_rotmat = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])

# The time range of interest
# look at a plot of dat.u versus dat.mpltime and can visually see where the data should be cropped
t_range = [
    # The instrument was in place starting at 12:08:30 on June 12,
    # 2012.
    735729.448,
    # dt.date2num(dt.datetime.datetime(2012, 6, 12, 12, 8, 30)),
    # The data is good to the end of the file.
    735731.305
]
# create a function for this that can determine the t_range using variance

# This is the filter to use for motion correction, defined
accel_filter = 0.1

# End user input section.
###############################

def t_range():
    # Read a file containing the adv data specified above:
    dat_raw = avm.read_nortek(fname)