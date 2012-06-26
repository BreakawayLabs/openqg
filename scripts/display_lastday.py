import scipy.io.netcdf as nc
import sys
import os

def convert_t_points(data):
    delta = data[1] - data[0]
    data = data - delta/2
    data = concatenate((data,  [data[-1] + delta]))
    return data

def plot_temp(title_str, temp_str, ncfile, max_x, max_y):
    title("%s mixed layer temperature ($\Delta$T (K))" % title_str)
    pcolor(convert_t_points(ncfile.variables["x"].data), convert_t_points(ncfile.variables["y"].data), ncfile.variables[temp_str].data, edgecolors='none', shading='flat', antialiased=False)
    clim(-35,35)
    colorbar()
    xlim(0, max_x)
    ylim(0, max_y)
    xlabel("E/W (km)")
    ylabel("N/S (km)")

def mixed_layer(aml, oml):

    max_x = max(aml.variables["x"].data)
    max_y = max(aml.variables["y"].data)

    clf()
    subplots_adjust(hspace=0.4)
    subplot(211)
    plot_temp("Atmosphere", "ast", aml, max_x, max_y)
    subplot(212)
    plot_temp("Ocean", "sst", oml, max_x, max_y)
    savefig("mixed_layer_temps.png")


def plot_pressure(title_str, k, ncfile, max_x, max_y, max_p):
    title("%s pressure" % title_str)
    pcolor(convert_t_points(ncfile.variables["x"].data), convert_t_points(ncfile.variables["y"].data), ncfile.variables["p"].data[k,:,:], edgecolors='none', shading='flat', antialiased=False)
    clim(-max_p,max_p)
    colorbar()
    xlim(0, max_x)
    ylim(0, max_y)
    xlabel("E/W (km)")
    ylabel("N/S (km)")

def pressure(atm, ocn):

    max_x = max(atm.variables["x"].data)
    max_y = max(atm.variables["y"].data)

    clf()
    figure(figsize=(8,12))
    subplots_adjust(hspace=0.6)
    subplot(611)
    plot_pressure("Atmosphere", 2, atm, max_x, max_y, 10000)
    subplot(612)
    plot_pressure("Atmosphere", 1, atm, max_x, max_y, 10000)
    subplot(613)
    plot_pressure("Atmosphere", 0, atm, max_x, max_y, 10000)
    subplot(614)
    plot_pressure("Ocean", 0, ocn, max_x, max_y, 50)
    subplot(615)
    plot_pressure("Ocean", 1, ocn, max_x, max_y, 50)
    subplot(616)
    plot_pressure("Ocean", 2, ocn, max_x, max_y, 50)
    savefig("pressures.png")
    #show()

if __name__ == '__main__':
    outdir = sys.argv[1]
    for filename in os.listdir(outdir):
        if "atmos_ml" in filename:
            atmos_ml = filename
        if "ocean_ml" in filename:
            ocean_ml = filename
        if "atmos_qg" in filename:
            atmos_qg = filename
        if "ocean_qg" in filename:
            ocean_qg = filename

    try:
        aml = nc.netcdf_file(os.path.join(outdir, atmos_ml))
        oml = nc.netcdf_file(os.path.join(outdir, ocean_ml))
        atm = nc.netcdf_file(os.path.join(outdir, atmos_qg))
        ocn = nc.netcdf_file(os.path.join(outdir, ocean_qg))
    except:
        print "Could not find all lastday files for mixed layer display"
        sys.exit(0)

    import matplotlib as mpl
    mpl.use('Agg')
    from pylab import *

    mixed_layer(aml, oml)

    pressure(atm, ocn)
