import os
import scipy.io.netcdf as nc
import matplotlib as mpl
mpl.use('Agg')
from pylab import *

def single_plot(time, var):
    title("%s (%s)" % (var.long_name, var.units))
    legend()
    plot(time, var.data)    

def plot_qg_mon(qg, output_file, title_str):

    ncfile = nc.netcdf_file(qg)
    print sorted(ncfile.variables.keys())
    time = ncfile.variables["time"].data
    pavg = ncfile.variables["pavg"].data
    qavg = ncfile.variables["qavg"].data

    print "BLARG", ncfile.variables["pavg"].units
    print "BLARG", ncfile.variables["pavg"].long_name

    clf()
    figure(figsize=(16,12))
    subplot(421)
    single_plot(time, ncfile.variables["pavg"])
    subplot(422)
    single_plot(time, ncfile.variables["qavg"])
    subplot(423)
    single_plot(time, ncfile.variables["circ"])
    subplot(424)
    single_plot(time, ncfile.variables["ctot"])
    subplot(425)
    single_plot(time, ncfile.variables["sfmax"])
    subplot(426)
    single_plot(time, ncfile.variables["sfmin"])
    subplot(427)
    single_plot(time, ncfile.variables["val"])
    xlabel("Time (years)")
    subplot(428)
    single_plot(time, ncfile.variables["entm"])
    xlabel("Time (years)")

    savefig(output_file)

if __name__ == '__main__':
    outdir = sys.argv[1]
    qgo = None
    qga = None
    for filename in os.listdir(outdir):
        if "qgo_mon" in filename:
            qgo = filename
        if "qga_mon" in filename:
            qga = filename
    print "OUTDIR", outdir
    print "QGO", qgo
    print "QGA", qga
    if qgo:
        plot_qg_mon(os.path.join(outdir, qgo), "qgo_mon.png", "Ocean")
    if qga:
        plot_qg_mon(os.path.join(outdir, qga), "qga_mon.png", "Atmosphere")
