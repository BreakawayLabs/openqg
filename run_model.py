"""
Main script for running Q-GCM models. Run

  python run_model.py --help

for uasge.
"""

import argparse
import subprocess
import sys
import os
import shutil

### Constant configuration classes for example experiments

class Config:

    glam = None
    ocn = None
    atm = None
    wind = None
    clock = None
    data = None

class DoubleGyre(Config):
    
    glam = "glam-dg.cdl"
    ocn = "ocn_dg.cdl"
    atm = "atm_dg.cdl"
    wind = "windstress-1.3.cdl"
    clock = "clock-0-5yr-6.cdl"

class DoubleGyreFast(DoubleGyre):

    clock = "clock-0-10dy-6.cdl"

class SouthernOceanOnly(Config):

    glam = "glam-so-only.cdl"
    ocn = "ocn_so_only.cdl"
    atm = None
    wind = "windstress-avges.cdl"
    clock = "clock-0-5yr-3.cdl"
    data = ["soctopog.26deg.10km.new.nc", "avges.nc"]

class SouthernOceanOnlyFast(SouthernOceanOnly):

    clock = "clock-0-10dy-3.cdl"
    
class SouthernOcean(Config):

    glam = "glam-so.cdl"
    ocn = "ocn_so.cdl"
    atm = "atm_so.cdl"
    wind = "windstress-1.5.cdl"
    clock = "clock-100-5yr-3.cdl"
    data = ["soctopog.26deg.10km.new.nc", 
            "lastday_so_ocn.nc", 
            "lastday_so_atm.nc",
            "lastday_so_aml.nc",
            "lastday_so_oml.nc"]

class SouthernOceanTau(SouthernOcean):

    wind = "windstress-1.5-udiff.cdl"

class SouthernOceanFast(SouthernOcean):

    clock = "clock-100yr-10dy-3.cdl"

class SouthernOceanTauFast(SouthernOceanTau):

    clock = "clock-100yr-10dy-3.cdl"

CONFIGS = {"dg": DoubleGyre,
           "dg_fast": DoubleGyreFast,
           "so_only": SouthernOceanOnly,
           "so_only_fast": SouthernOceanOnlyFast,
           "so": SouthernOcean,
           "so_fast": SouthernOceanFast,
           "so_tau": SouthernOceanTau,
           "so_tau_fast": SouthernOceanTauFast,
           }

###

def _mkdir(dirname):
    """ Equiv: mkdir -p dirname. """
    try:
        os.mkdir(dirname)
    except OSError:
        pass

def setup_from_config(config, setup_dir):
    """
    Take a Config object and generate an input directory.
    """
    _mkdir(setup_dir)
    shutil.copy("data/glam/%s" % config.glam,  "%s/glam.cdl" % setup_dir)
    shutil.copy("data/basin/%s" % config.ocn, "%s/ocn_basin.cdl" % setup_dir)
    if config.atm:
        shutil.copy("data/basin/%s" % config.atm, "%s/atm_basin.cdl" % setup_dir)
    shutil.copy("data/windstress/%s" % config.wind, "%s/windstress.cdl" % setup_dir)
    shutil.copy("data/clocks/%s" % config.clock, "%s/clock.cdl" % setup_dir)
    if config.data:
        for filename in config.data:
            shutil.copy("data/nc_files/%s" % filename, setup_dir)

def setup_from_dir(input_dir, setup_dir):
    """ Copy all .nc and .cdl files form input_dir to setup_dir. """
    _mkdir(setup_dir)
    for filename in os.listdir(input_dir):
        if filename.endswith((".nc", ".cdl")):
            shutil.copy("%s/%s" % (input_dir, filename), setup_dir)

stage_from_input = setup_from_dir

def _ncgen(filename, stage_dir):
    """
    Takes a filename <foo.cdl> in stage_dir and generates <foo.nc> in stage_dir
    by calling out to the ncgen shell command.
    """
    base = ".".join(filename.split(".")[:-1])
    cmd = "ncgen %s/%s.cdl -o %s/%s.nc" % (stage_dir, base, stage_dir, base)
    print cmd
    os.system(cmd)

def ncgen_all(stage_dir):
    """ Run ncgen on all .cdl in the directory directory. """
    for filename in os.listdir(stage_dir):
        if filename.endswith(".cdl"):
            _ncgen(filename, stage_dir)

def run_model(stage_dir, out_dir, exe, num_threads):
    # Run the model, piping stdout from the model to stdout of this script
    _mkdir(out_dir)
    cmd = "OMP_NUM_THREADS=%d %s %s %s" % (num_threads, exe, stage_dir, out_dir)
    #cmd = "valgrind %s %s %s" % (exe, stage_dir, out_dir)
    print cmd
    f = open("%s/stdout.txt" % out_dir, "w")
    p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p1.stdout
    while True:
        line = out.readline()
        if line:
            print line,
            f.write(line)
        else:
            break
    f.close()
    ret = p1.wait()

    if ret != 0:
        print "FAILED in run. Failed stage output:", out_dir
        sys.exit(ret)

def run_single(input_dir, stage_dir, out_dir, exe, num_threads):

    # cp input_dir/*.{nc,cdl} stage_dir
    stage_from_input(input_dir, stage_dir)

    # ncgen stage_dir/*cdl -> stage_dir/*nc
    ncgen_all(stage_dir)

    # q-gcm stage_dir -> out_dir
    run_model(stage_dir, out_dir, exe, num_threads)


def update_clock(filename):
    """
    Edit a clock file, updating 'tini' by adding 'trun'.
    """
    if not os.path.exists(filename):
        ncfile = filename.rsplit(".", 1)[0] + ".nc"
        subprocess.call("ncdump %s > %s" % (ncfile, filename), shell=True)
    for line in open(filename):
        if line.strip().startswith("tini ="):
            tini = float(line.strip().split('=')[1][:-1])
        if line.strip().startswith("trun ="):
            trun = float(line.strip().split('=')[1][:-1])
    tini += trun
    f = open(filename + ".tmp", "w")
    for line in open(filename):
        if "tini" in line and "=" in line:
            line = "\ttini = %f;\n" % tini
        if "trun" in line and "=" in line:
            line = "\ttrun = %f;\n" % trun
        f.write(line)
    f.close()
    os.rename("%s.tmp" % filename, filename)

def update_basin(filename):
    """
    Take a basin .cdl file and rewrite the stage values to 
    look for the appropriate "*_lastday.nc" file.
    """
    if not os.path.exists(filename):
        base = ".".join(filename.split(".")[:-1])
        cmd = "ncdump %s.nc > %s.cdl" % (base, base)
        print cmd
        os.system(cmd)
    if not os.path.exists(filename):
        print "NCDUMP FAILED"
        return
    for line in open(filename):
        if line.strip().startswith("name ="):
            name = line.strip().split('"')[1].strip()
            break
    f = open(filename + ".tmp", "w")
    for line in open(filename):
        if "state" in line and "=" in line:
            new_state = "%s_%s_lastday.nc" % (name, line.split("_")[0].strip().replace("oml", "ml").replace("aml", "ml"))
            old = line.split('"')
            old[1] = new_state
            new_line = '"'.join(old)
            f.write(new_line)
        else:
            f.write(line)
    f.close()
    os.rename("%s.tmp" % filename, filename)

def create_next_input_dir(stage_dir, out_dir, next_in):
    # take stage-n, out-n and create in-(n+1)
    _mkdir(next_in)
    files_to_copy = ["%s/%s" % (stage_dir, filename) for filename in os.listdir(stage_dir)]
    for filename in os.listdir(out_dir):
        if filename.endswith("lastday.nc"):
            files_to_copy.append("%s/%s" % (out_dir, filename))
    for filename in files_to_copy:
        shutil.copy(filename, next_in)

    update_basin("%s/ocn_basin.cdl" % next_in)
    update_basin("%s/atm_basin.cdl" % next_in)
    update_clock("%s/clock.cdl" % next_in)

def main(exe, output_dir, num_threads, repeats, input_dir, exp):
    
    # Create the base output directory
    _mkdir(output_dir)

    # Setup the very first directory, in-orig
    setup_dir = "%s/in-orig" % output_dir
    if input_dir is not None:
        # cp input_dir/*.{nc,cdl} -> setup_dir
        setup_from_dir(input_dir, setup_dir)
    elif exp is not None:
        setup_from_config(CONFIGS[exp], setup_dir)
    else:
        assert False

    # Copy everything from the input director into in-0
    shutil.copytree("%s/in-orig" % output_dir, "%s/in-0" % output_dir)

    for run in range(repeats):
        in_dir = "%s/in-%d" % (output_dir, run)
        stage_dir = "%s/stage-%d" % (output_dir, run)
        out_dir = "%s/out-%d" % (output_dir, run)
        next_in = "%s/in-%d" % (output_dir, run+1)

        # in-n -> stage-n -> out-n
        run_single(in_dir, stage_dir, out_dir, exe, num_threads)

        # (stage-n, out-n) -> in-(n+1)
        create_next_input_dir(stage_dir, out_dir, next_in)


def _strip_slash(dirname):
    if dirname is not None:
        if dirname.endswith("/"):
            dirname = dirname[:-1]
    return dirname

def parse_args():

    parser = argparse.ArgumentParser(description='Run a Q-GCM experiment.',
                                     epilog="One of the -e or -i options must be provided.")
    parser.add_argument('-x', '--exe', dest="exe", type=str, default="src/q-gcm", required=True,
                        help='Q-GCM executable to use')
    parser.add_argument("-o", "--output", dest="output_dir", type=str, required=True,
                        help='Output directory')

    parser.add_argument('-n', '--num_threads', dest="num_threads", type=int, default=4,
                        help='The number of threads to use when running in parallel')
    parser.add_argument('-r', '--runs', dest="repeats", type=int, default=1,
                        help='The number of sequential runs to perform')

    parser.add_argument('-i', '--input', dest="input_dir", type=str,
                        help='Directory containing input files.')
    parser.add_argument('-e', '--exp', dest="exp", type=str, choices=CONFIGS,
                        help='A predefined experiment.')

    args = parser.parse_args()

    if args.exp is None and args.input_dir is None:
        parser.print_usage()
        print >> sys.stderr, "%s: error: One of the -e or -i options must be provided." % sys.argv[0]
        sys.exit(1)

    print "EXE", args.exe
    print "OUTPUT", args.output_dir
    print "NUM THREADS", args.num_threads
    print "REPEATS", args.repeats
    print "INPUT", args.input_dir
    print "EXP", args.exp

    return args

if __name__ == '__main__':
    args = parse_args()
    main(args.exe, _strip_slash(args.output_dir), args.num_threads, args.repeats, _strip_slash(args.input_dir), args.exp)
