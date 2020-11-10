#!/usr/bin/python

## hydroplot.py
# Compare several different smoothing functions for the hydrophobicity of bR
# The bR structure is known so highlight the secondary structure features; from
#  http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=66360541

import argparse
import fasta_reader
import savitzky_golay

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
	#print('no display found. Using non-interactive Agg backend')
	mpl.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

def get_per_residue_values(sequence, table):
    values = []
    for residue in sequence:
        values.append(table[residue])
    return values

# Smooth a list of values using on a sliding filter.  The filter weights
# are in "weights" which must have an odd number of elements.  The middle
# element is the weight for the center of the window.
def smooth_values(values, weights):
    window_size = len(weights)
    if window_size%2 != 1:
        raise TypeError("smoothing requires an odd number of weights")
    
    half_window = (window_size-1)/2

    # Precompute the offset values for better performance.
    offsets = range(-half_window, half_window+1)
    offset_data = zip(offsets, weights)

    # normalize the weights in case the sum != 1
    total_weight = sum(weights)

    weighted_values = []
    
    for i in range(half_window, len(values)-half_window):
        weighted_value = 0.0
        for offset, weight in offset_data:
            weighted_value += weight*values[i+offset]
        weighted_values.append(weighted_value / total_weight)
    
    return weighted_values


def main(window, scale, infile, outfile, coord, dpi):
    # hydrophobicity indexes
    scaledesc = {'kd': 'Kyte&Doolittle hydrophobicity', 'hw': 'Hopp&Woods hydrophilicity'}
    scaledata = {}
    scaledata['kd'] = 	{ 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
		          'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
		          'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
		          'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
    scaledata['hw'] = 	{ 'A':-0.5,'R': 3.0,'N': 0.2,'D': 3.0,'C':-1.0,
		          'Q': 0.2,'E': 3.0,'G': 0.0,'H':-0.5,'I':-1.8,
		          'L':-1.8,'K': 3.0,'M':-1.3,'F':-2.5,'P': 0.0,
		          'S': 0.3,'T':-0.4,'W':-3.4,'Y':-2.3,'V':-1.5 }

#https://www.meme.net.au/savitzky-golay.html  
#aa	aa	Kyte-Doolittle	Hopp-Woods	Cornette	Eisenberg	Rose	Janin	Engelman GES
#A	Alanine		1.80	-0.50	0.20	0.62	0.74	0.30	1.60
#C	Cysteine	2.50	-1.00	4.10	0.29	0.91	0.90	2.00
#D	Aspartic acid	-3.50	3.00	-3.10	-0.90	0.62	-0.60	-9.20
#E	Glutamic acid	-3.50	3.00	-1.80	-0.74	0.62	-0.70	-8.20
#F	Phenylalanine	2.80	-2.50	4.40	1.19	0.88	0.50	3.70
#G	Glycine		-0.40	0.00	0.00	0.48	0.72	0.30	1.00
#H	Histidine	-3.20	-0.50	0.50	-0.40	0.78	-0.10	-3.00
#I	Isoleucine	4.50	-1.80	4.80	1.38	0.88	0.70	3.10
#K	Lysine		-3.90	3.00	-3.10	-1.50	0.52	-1.80	-8.80
#L	Leucine		3.80	-1.80	5.70	1.06	0.85	0.50	2.80
#M	Methionine	1.90	-1.30	4.20	0.64	0.85	0.40	3.40
#N	Asparagine	-3.50	0.20	-0.50	-0.78	0.63	-0.50	-4.80
#P	Proline		-1.60	0.00	-2.20	0.12	0.64	-0.30	-0.20
#Q	Glutamine	-3.50	0.20	-2.80	-0.85	0.62	-0.70	-4.10
#R	Arginine	-4.50	3.00	1.40	-2.53	0.64	-1.40	-12.3
#S	Serine		-0.80	0.30	-0.50	-0.18	0.66	-0.10	0.60
#T	Threonine	-0.70	-0.40	-1.90	-0.05	0.70	-0.20	1.20
#V	Valine		4.20	-1.50	4.70	1.08	0.86	0.60	2.60
#W	Tryptophan	-0.90	-3.40	1.00	0.81	0.85	0.30	1.90
#Y	Tyrosine	-1.30	-2.30	3.20	0.26	0.76	-0.40	-0.70
 

    average_weights = [1]*window
    triangle_weights = []
    for i in range(1, (((window-1)/2)+2)):
       triangle_weights.append(i)
    for i in reversed(range(1, (((window-1)/2)+1))):
       triangle_weights.append(i)

    # Savitzky-Golay smooth quadratic weights
    savitzky_golay_weights = [
        -60.15037594, -22.55639098, 10.61477222, 39.36311367,
        63.68863335, 83.59133127, 99.07120743, 110.12826183,
        116.76249447, 118.97390535, 116.76249447, 110.12826183,
        99.07120743, 83.59133127, 63.68863335, 39.36311367,
        10.61477222, -22.55639098, -60.15037594]

    savitzky_golay_weights = savitzky_golay.savitzky_golay_weights(window)
    
    record = fasta_reader.read_fasta_record(infile)
    
    values = get_per_residue_values(record.sequence, scaledata[scale])

    average_values = smooth_values(values, average_weights)
    triangle_values = smooth_values(values, triangle_weights)
    savitzky_golay_values = smooth_values(values, savitzky_golay_weights)
    
    ## Make the plot
    
    half_window = (len(average_weights)-1)/2
    x_data = range(half_window, half_window+len(average_values))
    
    plot(x_data, average_values, linewidth=1.0, label="average")
    plot(x_data, triangle_values, linewidth=1.0, label="triangle")
    plot(x_data, savitzky_golay_values, linewidth=1.0, label="Savitzky-Golay")

    legend(prop={'size': 10})

    # Draw a reasonable cutoff for membrane prediction
    # Value of 1.6 taken from
    #   http://arbl.cvmbs.colostate.edu/molkit/hydropathy/
#    axhline(y=1.6,color="gray")
    axhline(y=0,color="gray")
    
    # Draw the background mark
    if coord:
	    for clist in coord:
    		if (clist[0] > 0) and (clist[1] > clist[0]):
			axvspan(clist[0], clist[1], facecolor="yellow", alpha=0.4)
    
    # Show exactly the length of the sequence
    axis(xmin = 1, xmax = len(record.sequence))
    
    xlabel("residue number")
    ylabel(scaledesc[scale])
    
    record_id = record.title.split(" ", 1)[0]
    title("Protein " + record_id)
    
    if outfile:
    	savefig(outfile, dpi=dpi)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( '-w', '--window', type=int, choices=range(5,19,2), help='Window size',default=19)
    parser.add_argument( '-s', '--scale', type=str, choices=['kd','hw'], help='Available scales: kd (Kyte&Doolittle); hw (Hopp&Woods)',default='hw')
    parser.add_argument( '-i', '--infile', required=True, type=argparse.FileType('r'), help='Input protein sequence file',default=sys.stdin, nargs='?')
    parser.add_argument( '-o', '--outfile', required=True, type=argparse.FileType('w'), help="Output image plot file", default=sys.stdout, nargs='?')
    parser.add_argument( '-d', '--dpi', type=int, choices=range(80,301,1), help='DPI',default=80)
    parser.add_argument( '-c', '--coord', nargs=2, action='append', required=False, type=int, help="start,end coordinates for range mark")
    
    args = parser.parse_args()


#    parser.print_help()
    main(args.window,args.scale,args.infile,args.outfile,args.coord,args.dpi)
