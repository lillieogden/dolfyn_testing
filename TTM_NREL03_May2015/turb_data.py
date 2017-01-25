# To get started first import the DOLfYN ADV advanced programming
# interface (API):
import dolfyn.adv.api as avm
import dolfyn.adv.turbulence as turb


# Import matplotlib tools for plotting the data:
from matplotlib import pyplot as plt
import matplotlib.dates as dt
import numpy as np

##############################
# User input and customization

# The file to load:
fname = '/Users/lillie/turbulence_data/raw_data/TTM_NREL03_May2015.VEC'
# This file is available at:
# http://goo.gl/yckXtGk

# This is the vector from the ADV head to the body frame, in meters,
# in the ADV coordinate system.
body2head_vec = np.array([9.75, 2, -5.75]) * 0.0254


# This is the orientation matrix of the ADV head relative to the body.
# In this case the head was aligned with the body, so it is the
# identity matrix:
body2head_rotmat = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])

# The time range of interest.
t_range = [
    # The instrument was in place starting at 12:08:30 on June 12,
    # 2012.
    dt.date2num(dt.datetime.datetime(2012, 6, 12, 12, 8, 30)),
    # The data is good to the end of the file.
    np.inf
]

# This is the filter to use for motion correction:
accel_filter = 0.1

# End user input section.
###############################

# Read a file containing adv data:
dat_raw = avm.read_nortek(fname)

# Crop the data for t_range using DOLfYN's 'subset' method (creates a
# copy):
t_range_inds = (t_range[0] < dat_raw.mpltime) & (dat_raw.mpltime < t_range[1])
dat = dat_raw.subset(t_range_inds)
dat.props['body2head_vec'] = body2head_vec
dat.props['body2head_rotmat'] = body2head_rotmat

# Then clean the file using the Goring+Nikora method:
avm.clean.GN2002(dat)

# create the turbulence data object
turb_dat = turb.calc_turbulence(dat, 10)

# create the figure to add the turbulence data to in order to create a spectra
fig = plt.figure(1, figsize=[8, 4])
fig.clf()
ax = fig.add_axes([.14, .14, .8, .74])

# plot the turbulence data
ax.plot(turb_dat.mpltime, turb_dat.u, 'g-', rasterized=True)
bads = np.abs(turb_dat.u - dat_raw.u[t_range_inds])