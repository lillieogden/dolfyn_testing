# To get started first import the DOLfYN ADV advanced programming
# interface (API):
import dolfyn.adv.api as avm
from crop_data import t_set


# Import matplotlib tools for plotting the data:
from matplotlib import pyplot as plt
import matplotlib.dates as dt
import numpy as np


def crop_clean(t_range, dat_raw, body2head_rotmat, body2head_vec):
    # Crop the data for t_range using DOLfYN's 'subset' method (creates a
    # copy):
    t_range_inds = (t_range[0] < dat_raw.mpltime) & (dat_raw.mpltime < t_range[1])
    dat = dat_raw.subset(t_range_inds)
    dat.props['body2head_vec'] = body2head_vec
    dat.props['body2head_rotmat'] = body2head_rotmat

    # Then clean the file using the Goring+Nikora method:
    avm.clean.GN2002(dat)
    dat_cln = dat.copy()
    return dat, dat_cln, t_range_inds


def motion(dat, dat_cln, accel_filter):
    """Perform motion correction (including rotation into earth frame):"""

    avm.motion.correct_motion(dat, accel_filter)

    # Rotate the uncorrected data into the earth frame,
    # for comparison to motion correction:
    avm.rotate.inst2earth(dat_cln)

    # ax.plot(dat.mpltime, dat.u, 'b-')

    # Then rotate it into a 'principal axes frame':
    avm.rotate.earth2principal(dat)
    avm.rotate.earth2principal(dat_cln)

    # save the data
    dat.save('TTM_NREL03_May2015_rotated2principal.h5')

# And reload the data:
# dat_bin_copy = avm.load('TTM_NREL03_May2015_rotated2principal.h5')


def crop_plot(dat_raw, dat, t_range_inds, t_range):
    """Create a figure for comparing screened data to the original."""

    fig = plt.figure(1, figsize=[8, 4])
    fig.clf()
    ax = fig.add_axes([.14, .14, .8, .74])

    # Plot the raw (unscreened) data:
    # this data has not been cropped or cleaned
    # first, convert the num_time to date_time, and plot this versus dat_raw.u
    date_time_raw = dt.num2date(dat_raw.mpltime)
    ax.plot(date_time_raw, dat_raw.u, 'r-', rasterized=True)

    # Plot the screened data:
    date_time = dt.num2date(dat.mpltime)
    ax.plot(date_time, dat.u, 'g-', rasterized=True)
    bads = np.abs(dat.u - dat_raw.u[t_range_inds])
    ax.text(0.55, 0.95,
            "%0.2f%% of the data were 'cleaned'\nby the Goring and Nikora method."
            % (np.float(sum(bads > 0)) / len(bads) * 100),
            transform=ax.transAxes,
            va='top',
            ha='left',
            )

    # Add some annotations:
    ax.axvspan(dat_raw.mpltime[0], t_range[0], zorder=-10, facecolor='0.9', edgecolor='none')
    ax.text(0.13, 1.0, 'Mooring falling\ntoward seafloor', ha='center', va='top', transform=ax.transAxes, size='small')
    ax.text(0.3, 0.6, 'Mooring on seafloor', ha='center', va='top', transform=ax.transAxes, size='small')
    ax.annotate('', (0.25, 0.4), (0.4, 0.4), arrowprops=dict(facecolor='black'))

    # Finalize the figure
    # Format the time axis:
    # tkr = dt.MinuteLocator(interval=5)
    # frmt = dt.DateFormatter('%H:%M')
    # ax.xaxis.set_major_locator(tkr)
    # ax.xaxis.set_minor_locator(dt.MinuteLocator(interval=1))
    # ax.xaxis.set_major_formatter(frmt)
    # ax.set_ylim([-3, 3])

    # Label the axes:
    ax.set_ylabel('$u\,\mathrm{[m/s]}$', size='large')
    ax.set_xlabel('Time [June 12, 2012]')
    ax.set_title('Data cropping and cleaning')
    # ax.set_xlim([dt.date2num(dt.datetime.datetime(2012, 6, 12, 12)),
    #              dt.date2num(dt.datetime.datetime(2012, 6, 12, 12, 30))])

    # Save the figure:
    fig.savefig('/Users/lillie/turbulence_data/plots/TTM_NREL03_May2015/crop_data.pdf')
    # end cropping figure


def spectra(dat_bin, dat_cln_bin):
    """Create a plot for looking at the turbulence spectra"""
    # Figure to look at spectra
    fig2 = plt.figure(2, figsize=[6, 6])
    fig2.clf()
    ax = fig2.add_axes([.14, .14, .8, .74])

    ax.loglog(dat_bin.freq, dat_bin.Suu_hz.mean(0),
              'b-', label='motion corrected')
    ax.loglog(dat_cln_bin.freq, dat_cln_bin.Suu_hz.mean(0),
              'r-', label='no motion correction')

    # Add some annotations

    # adds the solid horizontal black line
    ax.axhline(1.7e-4, color='k', zorder=21)

    ax.text(2e-3, 1.7e-4, 'Doppler noise level', va='bottom', ha='left',)

    ax.text(1, 2e-2, 'Motion\nCorrection')
    ax.annotate('', (3.6e-1, 3e-3), (1, 2e-2),
                arrowprops={'arrowstyle': 'fancy',
                            'connectionstyle': 'arc3,rad=0.2',
                            'facecolor': '0.8',
                            'edgecolor': '0.6',
                            },
                ha='center',
                )

    ax.annotate('', (1.6e-1, 7e-3), (1, 2e-2),
                arrowprops={'arrowstyle': 'fancy',
                            'connectionstyle': 'arc3,rad=0.2',
                            'facecolor': '0.8',
                            'edgecolor': '0.6',
                            },
                ha='center',
                )

    # Finalize the figure
    ax.set_xlim([1e-3, 20])
    ax.set_ylim([1e-4, 1])
    ax.set_xlabel('frequency [hz]')
    ax.set_ylabel('$\mathrm{[m^2s^{-2}/hz]}$', size='large')

    f_tmp = np.logspace(-3, 1)
    ax.plot(f_tmp, 4e-5 * f_tmp ** (-5. / 3), 'k--')

    ax.set_title('Velocity Spectra')
    ax.legend()
    ax.axvspan(1, 16, 0, .2, facecolor='0.8', zorder=-10, edgecolor='none')
    ax.text(4, 4e-4, 'Doppler noise', va='bottom', ha='center',
            # bbox=dict(facecolor='w', alpha=0.9, edgecolor='none'),
            zorder=20)

    fig2.savefig('/Users/lillie/turbulence_data/plots/TTM_NREL03_May2015/motion_vel_spec.pdf')


def turb_energy(dat_cln_bin):
    fig = plt.figure(1, figsize=[8, 4])
    fig.clf()
    ax = fig.add_axes([.14, .14, .8, .74])

    # first, convert the num_time to date_time, and plot this versus dat_raw.u
    date_time = dt.num2date(dat_cln_bin.mpltime)

    # plot the data
    ax.plot(date_time, dat_cln_bin.upup_, 'r-', rasterized=True)
    ax.plot(date_time, dat_cln_bin.vpvp_, 'g-', rasterized=True)
    ax.plot(date_time, dat_cln_bin.wpwp_, 'b-', rasterized=True)

    # label axes
    ax.set_xlabel('Time')
    ax.set_ylabel('Turbulent Energy $\mathrm{[m^2/s^2]}$', size='large')

    fig.savefig('/Users/lillie/turbulence_data/plots/TTM_NREL03_May2015/turb_energy_spec.pdf')


def reynolds_stress(dat_cln_bin):
    fig = plt.figure(1, figsize=[8, 4])
    fig.clf()
    ax = fig.add_axes([.14, .14, .8, .74])

    # first, convert the num_time to date_time, and plot this versus dat_raw.u
    date_time = dt.num2date(dat_cln_bin.mpltime)

    # plot the data
    ax.plot(date_time, dat_cln_bin.upvp_, 'r-', rasterized=True)
    ax.plot(date_time, dat_cln_bin.upwp_, 'g-', rasterized=True)
    ax.plot(date_time, dat_cln_bin.vpwp_, 'b-', rasterized=True)

    # label axes
    ax.set_xlabel('Time')
    ax.set_ylabel('Reynolds Stresses $\mathrm{[m^2/s^2]}$', size='large')

    fig.savefig('/Users/lillie/turbulence_data/plots/TTM_NREL03_May2015/reynolds_stress.pdf')


def main():
    #############################
    # User input and customization, global variables

    # The file to load:
    fname = '/Users/lillie/turbulence_data/raw_data/TTM_NREL03_May2015.VEC'
    name = 'TTM_NREL03_May2015'
    # Read a file containing the adv data specified above:
    dat_raw = avm.read_nortek(fname)

    # This is the vector from the ADV head to the body frame, in meters, in the ADV coordinate system.
    body2head_vec = np.array([9.75, 2, -5.75]) * 0.0254
    # = np.array([0.48,-0.07, -0.27])

    # This is the orientation matrix of the ADV head relative to the body.
    # In this case the head was aligned with the body, so it is the identity matrix:
    body2head_rotmat = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])
    # = np.eye(3)

    # the time range of interest
    # 1. look at a plot of dat.u versus dat.mpltime and can visually see where the data should be cropped
    # t_range = [735729.450, 735731.305]

    # 2. use variance thresholds to determine the t_range
    t_range = t_set(dat_raw, n_bin=8000, n_fft=4096)

    # This is the filter to use for motion correction, defined
    accel_filter = 0.1

    # End user input section.
    ###############################

    # crop and clean the data
    dat, dat_cln, t_range_inds = crop_clean(t_range, dat_raw, body2head_rotmat, body2head_vec)

    # perform motion correction
    motion(dat, dat_cln, accel_filter)

    # create a velocity distribution, illustrating the cropping and cleaning process
    crop_plot(dat_raw, dat, t_range_inds, t_range)

    # Average the data and compute turbulence statistics
    dat_bin = avm.calc_turbulence(dat, n_bin=19200,
                                  n_fft=4096)
    dat_cln_bin = avm.calc_turbulence(dat_cln, n_bin=19200,
                                      n_fft=4096)

    # make graphs of the velocity spectra, the tke and the reynolds stress
    spectra(dat_bin, dat_cln_bin)
    turb_energy(dat_cln_bin)
    reynolds_stress(dat_cln_bin)


# run the program
main()