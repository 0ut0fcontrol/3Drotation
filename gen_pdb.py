#!/usr/bin/env python2
import os,sys
from eft_calculator import EFT_calculator

calculator = EFT_calculator('wtr','alc')
#calculator = EFT_calculator('alc','wtr')
calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
root = calculator.com.frg + '_' + calculator.probe.frg + '.confs.dat'
if not os.path.exists(root):os.mkdir(root)
def mol2mol_init(ele):
    mol = [[i,0.0,0.0,0.0] for i in ele]
    return mol
size = 200
folder_id = 0
file_count = 0
#confs =  calculator.grid.gen_grid_x()
#for idx, coors in calculator.gen_PDB(confs): 
for idx, coors in calculator.gen_PDB(): 
#for id, coors in calculator.gen_atomic_coors(0,10): 
    #print(idx, coors)

    if  file_count%size == 0:
        folder = os.path.join(root,"EFT_%04d"%(folder_id))
        if not os.path.exists(folder):os.mkdir(folder)
        folder_id += 1
    pdb = open("%s/eft.%s.pdb"%(folder,str(idx)),"w")
    pdb.write(coors)
    pdb.close()
    file_count += 1

