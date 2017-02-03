#!/bin/python
import argparse
from math import cos, pi, sqrt

ARCSECOND_TO_DEG = 1./(60*60)
DEG_TO_RAD = 180. / pi

n_filters = 10


# FRED files have the following format:
## STANDARD MAGNITUDES ONLY
## FILCON ver 3.0
## RA (J2000)    DEC    CCDX      CCDY     Flags      HJD   Airmass  Set     Group    Field      Filt   Mag     Error    dmag    sys night
## 102.5140180   0.2598700  1828.670    16.950 0 0 56295.536090 2.133    1         92 0020110040     8  10.7481  0.0050  0.0460    31 56295

def read_fred_line(line):
    """Reads a line from a FRED file and outputs a properly-typed Python dict."""

    # Ideally we would use astropy.asciitable or numpy.loadtxt, but the
    # merged FRED file will be so big we cannot read it all at once. Instead
    # we tokenize and store the results in a Python dictionary.

    output = {}

    if line[0] == '#':
        return output

    line = line.strip()
    line = line.split()

    output['ra']     = float(line[0] )
    output['dec']    = float(line[1] )
    output['ccdx']   = float(line[2] )
    output['ccdy']   = float(line[3] )
    output['flag1']  = int(  line[4], 10)
    output['flag2']  = int(  line[5], 10)
    output['hjd']    = float(line[6] )
    output['avexx']  = float(line[7] )
    output['kset']   = int(  line[8], 10)
    output['group']  = int(  line[9], 10)
    output['star']   = int(  line[10], 10)
    output['filter_id'] = int(  line[11], 10) - 1 # filter indexing begins at 1
    output['xmag1']  = float(line[12])
    output['xerr1']  = float(line[13])
    output['dmag']   = float(line[14])
    output['sys']    = int(  line[15], 10)
    output['night']  = int(  line[16], 10)

    return output

def compute_statistics(tokenized_lines):
    """Computes the average RA, DEC, magnitude, and tracks the number of
    measurements of each filter for a list of tokenized FRED entries."""
    output = {}

    ra = 0
    ra2 = 0
    dec = 0
    dec2 = 0
    filter_id = 0
    filter_count = [0] * n_filters
    mags = [0] * n_filters
    mags2 = [0] * n_filters
    star_name = ''

    n_lines = len(tokenized_lines)
    for line in tokenized_lines:
        t_flag1 = line['flag1']
        t_xmag = line['xmag1']

        # skip measurements whose flags are set or magnitudes are
        # nonsensical
        if t_flag1 != 0 or t_xmag >= 90:
            continue

        # update the name, we needn't check that it matches
        star_name = line['star']

        # add in the RA and DEC
        t_ra = line['ra']
        t_dec = line['dec']
        ra += t_ra
        ra2 += t_ra * t_ra
        dec += t_dec
        ra2 += t_dec * t_dec

        # add in the magnitude
        filter_id = line['filter_id']
        filter_count[filter_id] += 1
        mags[filter_id] += t_xmag
        mags2[filter_id] += t_xmag * t_xmag

    # compute averages
    ra /= n_lines
    dec /= n_lines
#    raerr = sqrt((n_lines * ra2 - ra*ra) / (n_lines * n_lines - 1))
#    decerr = sqrt((n_lines * dec2 - dec*dec) / (n_lines * n_lines - 1))
#    for filter in range(0, n_filters):
#        fc = filter_count[filter]
#        if fc > 0:
#            mags[filter] /= fc

    # todo: return as a dictionary.
    output['star'] = star_name
    output['ra'] = ra
    output['ra2'] = ra2
    output['dec'] = dec
    output['dec2'] = dec2
    output['mags'] = mags
    output['filter_count'] = filter_count
    output['n_entries'] = n_lines

    return output

def write_statistics(stats, outfile):
#          write (3,9015) sname(k),ra(k),ra2(k),dec(k),dec2(k),
#     $     xnr(k),xmr(k),(smag(i,k),smag2(i,k),
#     $     serr(i,k),xnm(i,k),i=1,MAXFILT)
# 9015     format (a10,4d24.16,2i8,11(3e16.8,i10))
    star_name = stats['star']
    ra = stats['ra']
    ra2 = stats['ra2']
    dec = stats['dec']
    dec2 = stats['dec2']
    # TODO: figure these out:
    xnr = int(0)
    xmr = int(0)
    stats['errors'] = n_filters * [0]

    outline  = '{:0>10}'.format(star_name)
    outline += '{:24.16e}{:24.16e}{:24.16e}{:24.16e}'.format(ra, ra2, dec, dec2)
    outline += '{:>8}{:>8}'.format(xnr, xmr)
    for filter_id in range(0, n_filters):
        smag = stats['mags'][filter_id]
        smag2 = smag * smag
        serr = stats['errors'][filter_id]
        xnm = 0
        outline += '{:16.8e}{:16.8e}{:16.8e}{:>10}'.format(smag, smag2, serr, xnm)

    outfile.write(outline + "\n")

def combine_records(out_filename, in_filename, match_distance):
    """Parses the (sorted) input file to generate global statistics. Results
    are written in outfile."""

    ra_prev = 0
    dec_prev = 0
    entries = []

    outfile = open(out_filename, 'w')
    firstline = True;
    results = 0

    # read the file and parse the results
    with open(in_filename) as infile:
        for line in infile:

            # parse the line, skip lines with no values (comment lines)
            line = read_fred_line(line)
            if not line: # empty Python dicts return false
                continue

            # pull out the RA and DEC
            ra  = line['ra']
            dec = line['dec']

            if firstline:
                ra_prev = ra
                dec_prev = dec
                firstline = False

            # merge entries that are within the matching radius
            distance = ((ra - ra_prev)*cos(dec * DEG_TO_RAD))**2 + (dec - dec_prev)**2
            if distance < match_distance:
                entries.append(line)
            else:
                # We don't have match, compute average statistics
                stats = compute_statistics(entries)
                write_statistics(stats, outfile)

                # clear out the entries buffer and advance to the next star
                entries = []
                entries.append(line)

            # update the coordinates to match the latest star
            # (a hacky work-around for incorrectly split matches)
            ra_prev = ra
            dec_prev = dec


def sort_and_merge(tempdir, output, inputs, max_processes=8):
    """Sorts and merges the contents of `inputs` using gnuSort.
    Files are sorted by the combined key (ra,dec) (fields 1,2)"""
    import ntpath
    import subprocess
    import os

    processes = set()
    command = '/usr/bin/sort'
    sorted_files = []

    print("Start: Merging and Sorting Process")

    # Sort each file writing the results to the tempdir
    # executing gnuSort as `sort -n -o output input`
    # This is done in parallel using subprocess.Popen to leverage multi-core processors
    for filename in inputs:
        basename = ntpath.basename(filename)
        infile = filename
        outfile = tempdir + '/' + basename
        sorted_files.append(outfile)

        # notify the user
        print("Sorting " + infile)

        # run sort, following the mandated limits above
        processes.add(subprocess.Popen([command, '-n', '-k' '2,2', '-k', '1,1', '-o', outfile, infile]))
        # notify the user of progress
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

    # wait for all children processes to exit
    for p in processes:
        if p.poll() is None:
            p.wait()

    # merge the sorted files into a single output
    # i.e. `sort -m output inputs`
    print("Merging FRED files")
    cmd = [command, '-n', '-k', '2,2', '-k', '1,1', '-m', '-o', output]
    cmd.extend(sorted_files)
    subprocess.call(cmd)

    # remove the intermediate files
    print("Removing intermediate files")
 #   for filename in sorted_files:
 #       subprocess.call(['/usr/bin/rm', filename])

    print("Finish: Merging and Sorting")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('output', help="Output file containing the merged results")
    parser.add_argument('input', nargs='+', help="Input files which will be merged")
    parser.add_argument('--max-processes', type=int, default=8,
                        help="Maximum number of processes to use when sorting")
    parser.add_argument('--skip-sort', action='store_true',
                        help="Assume sorted.fred exists and skip the FRED sort/merge step")
    parser.add_argument('--radius', type=float, default=2.0,
                        help="Radius for matching stars in arseconds.")

    args = parser.parse_args()

    outdir = '/home/kloppenb/apass-source/scratch/'
    output = args.output
    inputs = args.input
    temp_dir = '/home/kloppenb/tmp/'
    sorted_fred = outdir + '/sorted.fred'

    # Stage 1: Sort individual FRED files and them merge them into a single,
    # sorted FRED file.
    if not args.skip_sort:
        sort_and_merge(temp_dir, sorted_fred, inputs, max_processes=args.max_processes)

    # Stage 2: Read in the sorted fred file and merge adjacent records
    match_radius = args.radius * ARCSECOND_TO_DEG
    match_distance = match_radius**2
    combine_records(output, sorted_fred, match_distance)

