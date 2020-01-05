import MDAnalysis as mda
import numpy as np
import sys
from MDAnalysis.analysis import align

top_file = sys.argv[1]


x_min = -40
x_max = 40
x_width = 0.1
num_x_bins = int((x_max-x_min)/x_width)
xcutoff = 40

y_min = -40
y_max = 40
y_width = 0.1
num_y_bins = int((y_max-y_min)/y_width)
ycutoff = 40

z_min = -10
z_max = 10
z_width = 0.2
num_z_bins = int((z_max-z_min)/z_width)
zcutoff = 10


hist = np.zeros((num_x_bins, num_y_bins, num_z_bins), dtype=float)

count = 0


i = 9
# loop over PDI molecules to center on 
while i <= 9:
    print "Centering on PDI %s" %(i)
    index = 0
    res_num = i
    # load in trajectory centered around the correct PDI 
    traj_file = "wrap_traj/RADA_stacked.run29.%s.dcd" %(res_num) 
    u = mda.Universe(top_file, traj_file)

    
    # Make selection of PDI core and amino acids. 
    # Atoms to exclude from PDI (glycine)
    ylg_atoms = ["HA1", "HA2", "CA1", "CA", "OA"]
    gly_atoms = ["HB1", "HB2", "CB1", "CB", "OB"]
    PDI_no_gly = u.select_atoms("resid %s and not (name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s)" %(res_num, ylg_atoms[0], ylg_atoms[1], ylg_atoms[2], ylg_atoms[3], ylg_atoms[4], gly_atoms[0], gly_atoms[1], gly_atoms[2], gly_atoms[3], gly_atoms[4]))
    wat_sel = u.select_atoms("resname WAT and name O")
    
    # Translate and align to the principal axes to create reference structure 
    u.atoms.translate(-PDI_no_gly.center_of_mass())
    u.atoms.rotate(PDI_no_gly.principal_axes())        

    # Write out a PDB of the reference structure
    with mda.Writer('peptide_pos/%s.align.pdb' %(i),u.atoms.n_atoms) as PDB:    
        PDB.write(u.atoms)


    # open an out file for the PDB containing residue center of mass locations
    out = open("peptide_pos/system.%s.pdb" %(res_num), 'w')


    # open reference structure in another pdb
    u2 = mda.Universe('peptide_pos/%s.align.pdb' %(i))
    ref = u2.select_atoms("resid %s and not (name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s or name %s)" %(res_num, ylg_atoms[0], ylg_atoms[1], ylg_atoms[2], ylg_atoms[3], ylg_atoms[4], gly_atoms[0], gly_atoms[1], gly_atoms[2], gly_atoms[3], gly_atoms[4]))
    step = 1
    
    with mda.Writer('peptide_pos/%s.rotated.dcd' %(i),u.atoms.n_atoms) as W:
        # loop over trajectory:
        for ts in u.trajectory:
            print step
            # align each frame to the reference structure
            u.atoms.translate(-PDI_no_gly.center_of_mass())
            R, rmsd = align.rotation_matrix(PDI_no_gly.positions, ref.positions)
            u.atoms.rotate(R)
            # write a dcd of the aligned structure
            W.write(u.atoms)
            # Store oxygen coordinates in histogram
            for atom in wat_sel.atoms:
                if np.abs(atom.position[0]) < xcutoff:
                    if np.abs(atom.position[1]) < ycutoff:
                        if np.abs(atom.position[2]) < zcutoff:

                            hist_x_bin = int((atom.position[0] - x_min)/x_width)
                            hist_y_bin = int((atom.position[1] - y_min)/y_width)
                            hist_z_bin = int((atom.position[2] - z_min)/z_width)
                #print atom.position, hist_x_bin,hist_y_bin,hist_z_bin
                #if 0 <= hist_x_bin <=num_x_bins-1 and  0 <= hist_y_bin <= num_y_bins-1 and 0 <= hist_z_bin <= num_z_bins-1:  
                # print "THIS ONE \n\n\n\n\n\n\n\n\n"
                            hist[hist_x_bin,hist_y_bin,hist_z_bin] += 1
                            count += 1
            step += 1          
    i += 17


out = open("RADA8_wat_gofr_3d_PDI9.dat",'w')
for x in range(num_x_bins):
    for y in range(num_y_bins):
        for z in range(num_z_bins):
            if hist[x,y,z] > 0: 
                out.write("  %5s  %5s  %5s  %10.5f  %s\n" %(x,y,z,hist[x,y,z],count))
