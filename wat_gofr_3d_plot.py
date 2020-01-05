import MDAnalysis as mda
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl


def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r

def plot_ram(data_combined,extent,system):
    #heatmap, xedges, yedges = np.histogram2d(psi, phi, bins=(nbins,nbins), range=[[-180, 180], [-180, 180]], normed=norm_flag)
    #heatmap *= norm_const
    # Plot heatmap
    plt.clf()
    #plt.axhline(0, color='k')
    #plt.axvline(0, color='k')
    plt.ylabel(r'')
    plt.xlabel(r'')
    my_cmap = plt.cm.get_cmap('seismic')
    my_cmapr = reverse_colourmap(my_cmap)
    my_cmapr.set_under('w')

    plt.imshow(data_combined, extent=extent, origin='lower', cmap = my_cmapr)
    ax, _ = mpl.colorbar.make_axes(plt.gca(), shrink=1.0)
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=my_cmapr,norm=mpl.colors.Normalize(vmin=0, vmax=5))
    cbar.set_clim(-4.0, 6.0)
    plt.clim(-4,6)
    plt.savefig('%s.wat_gofr.png' %(system))
    plt.savefig('%s.wat_gofr.eps' %(system))
    plt.savefig('%s.wat_gofr.pdf' %(system))



data = np.loadtxt(sys.argv[1])
out_file = sys.argv[2]

x_min = -40
x_max = 40
x_width = 0.1
num_x_bins = int((x_max-x_min)/x_width)

y_min = -40
y_max = 40
y_width = 0.1
num_y_bins = int((y_max-y_min)/y_width)

z_min = -10
z_max = 10
z_width = 0.2
num_z_bins = int((z_max-z_min)/z_width)

print " Done loading data"
wat_den = 0.0334 # atoms/A^3

z_cutoff = 1.5
z_cutoff_bin_max = int((z_cutoff-z_min)/z_width)
z_cutoff_bin_min = int((-z_cutoff-z_min)/z_width)
xy_cutoff = 30
xy_cutoff_bin_max = int((xy_cutoff-x_min)/x_width)
xy_cutoff_bin_min = int((-xy_cutoff-x_min)/x_width)

data_combined = np.zeros((num_x_bins,num_y_bins),dtype=float)


 
for i in range(len(data[:,0])):
    if z_cutoff_bin_min <= data[i,2] <= z_cutoff_bin_max:
#        if xy_cutoff_bin_min <= data[i,0] <= xy_cutoff_bin_max:
#            if xy_cutoff_bin_min <= data[i,1] <= xy_cutoff_bin_max:
        data_combined[int(data[i,0]),int(data[i,1])] += data[i,3]

print "About to plot"
Norm = wat_den*5000*30*x_width*z_width*z_cutoff
for i in range(num_x_bins):
    for j in range(num_y_bins):
        data_combined[i,j] /= Norm

extent = [x_min, x_max, y_min, y_max]

plot_ram(data_combined.T,extent,out_file)
