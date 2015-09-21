#!/usr/bin/python2.7.2
import argparse
import h5py
import numpy as np

#merge

parser = argparse.ArgumentParser(description='merge hdf5 files from perms2h5; this is NOT a general utility for merging h5 files')

parser.add_argument("-i","--in-files",type=str,nargs='*',help="Which files are to be merged")
parser.add_argument("-o","--out-file",type=str,help="Name of output file")
parser.add_argument("-n","--nperms",type=int,help="Number of perms read at one time")

args=parser.parse_args()

inf_list = args.in_files
outF = h5py.File(args.out_file,"w")
total_in = len(inf_list)

N=args.nperms

inF = h5py.File(inf_list[0],"r") 

for name in inF:
    outF.create_group(name)
    for grp in inF[name]:
        single_dim = inF[name+'/'+grp].shape
        if grp == 'permutations' : 
            total_dim = [a for a in single_dim]
            total_dim[0]=total_dim[0]*total_in
            total_dim = tuple(total_dim)
            outF[name].create_dataset(grp,total_dim,dtype=inF[name+'/'+grp].dtype,chunks=inF[name+'/'+grp].chunks,compression="gzip",compression_opts=9)
        else:
            outF[name].create_dataset(grp,single_dim,dtype=inF[name+'/'+grp].dtype,data=inF[name+'/'+grp].value,chunks=single_dim,compression="gzip",compression_opts=9)

a = 0 
for a in range(len(inf_list)):
    if a > 0:
        inF.close()
        inF= h5py.File(inf_list[a],"r")
        
    single_dim = inF['/Perms/permutations'].shape
    i,j= 0,N
    while j<= single_dim[0]:
        outF['/Perms/permutations'][(i+a*single_dim[0]):(j+a*single_dim[0])]= inF['/Perms/permutations'][i:j]
        outF.flush()
        i += N
        j += N

inF.close()
outF.close()




     
                               
        



