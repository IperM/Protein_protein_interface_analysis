#!/usr/bin/env python
# coding: utf-8

# # Structure checking tutorial
# 
# A complete checking analysis of a single structure follows.
# use .revert_changes() at any time to recover the original structure

# Structure checking is a key step before setting up a protein system for simulations. 
# A number of normal issues found in structures at Protein Data Bank may compromise the success of the simulation, or may suggest that longer equilibration procedures are necessary.
# 
# The biobb_structure_checking modules allow to 
# - Do basic manipulations on structures (selection of models, chains, alternative locations
# - Detect and fix amide assignments, wrong chirality
# - Detect and fix protein backbone issues (missing fragments, and atoms, capping)
# - Detect and fix missing side-chain atoms
# - Add hydrogen atoms according to several criteria
# - Detect and classify clashes
# - Detect possible SS bonds
# 
# biobb_structure_checking modules can used at the command line biobb_structure_checking/bin/check_structure
# 

# In[1]:

# ## Installation

# #### Basic imports and initialization

# In[2]:


import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
base_dir_path=biobb_structure_checking.__path__[0]
args = cts.set_defaults(base_dir_path,{'notebook':True})


# ## General help

# In[3]:


with open(args['commands_help_path']) as help_file:
    print(help_file.read())
#TODO: prepare a specific help method
# print_help(command)


# Set input (PDB or local file, pdb or mmCif formats allowed) and output (local file, pdb format).  
# Use pdb:pdbid for downloading structure from PDB (RCSB)

# In[4]:


base_path = '/home/pablo/Documentos/BIOPHYSICS_PROJECT/BioPhysics/Notebooks/ '
args['input_structure_path'] = base_path + '7ekf.cif'
args['output_structure_path'] = base_path + '7ekf_fixed.pdb'
args['output_structure_path_charges'] = base_path + '7ekf_fixed.pdbqt'
args['debug'] = False
args['verbose'] = False


# Initializing checking engine, loading structure and showing statistics

# In[5]:



st_c = StructureChecking(base_dir_path, args)


# #### models
# Checks for the presence of models in the structure. 
# MD simulations require a single structure, although some structures (e.g. biounits) may be defined as a series of models, in such case all of them are usually required.  
# Use models('--select N') to select model num N for further analysis

# In[6]:


st_c.models()


# #### chains
# Checks for chains (also obtained from print_stats), and allow to select one or more.   
# MD simulations are usually performed with complete structures. However input structure may contain several copies of the system, or contains additional chains like peptides or nucleic acids that may be removed. 
# Use chains('X,Y') to select chain(s) X and Y to proceed

# In[7]:


st_c.chains()


# #### altloc
# Checks for the presence of residues with alternative locations. Atoms with alternative coordinates and their occupancy are reported.  
# MD simulations requires a single position for each atom.  
# Use altloc('occupancy | alt_ids | list of res:id) to select the alternative
# 

# In[8]:


st_c.altloc()


# We need to choose one of the alternative forms for each residue

# In[9]:


st_c.altloc('occupancy')


# In[10]:


st_c.altloc()


# #### metals
# Detects HETATM being metal ions allow to selectively remove them.  
# To remove use metals (' All | None | metal_type list | residue list ')

# In[11]:


st_c.metals()


# #### ligands
# Detects HETATM (excluding Water molecules) to selectively remove them.  
# To remove use ligands('All | None | Residue List (by id, by num)')
# 

# In[12]:


st_c.ligands()


# In[13]:


st_c.ligands('All')


# In[14]:


st_c.ligands()


# #### rem_hydrogen
# Detects and remove hydrogen atoms. 
# MD setup can be done with the original H atoms, however to prevent from non standard labelling, remove them is safer.  
# To remove use rem_hydrogen('yes')
# 

# In[15]:


st_c.rem_hydrogen()


# #### water
# Detects water molecules and allows to remove them
# Crystallographic water molecules may be relevant for keeping the structure, however in most cases only some of them are required. These can be later added using other methods (titration) or manually.
# 
# To remove water molecules use water('yes')
# 

# In[16]:


st_c.water()


# In[17]:


st_c.water("yes")


# #### amide
# Amide terminal atoms in Asn ang Gln residues can be labelled incorrectly.  
# amide suggests possible fixes by checking the sourrounding environent.
# 
# To fix use amide ('All | None | residue_list')
# 
# Note that the inversion of amide atoms may trigger additional contacts. 

# In[18]:


st_c.amide()


# Fix all amide residues and recheck

# In[19]:


st_c.amide('all')


# Comparing both checks it becomes clear that GLN A42, GLN E498, ASN A103, and ASN A194 are getting new contacts as thay have both changed, ASN E394 is worse as it has now two contacts

# In[20]:


st_c.amide('A42,A103')


# In[21]:


st_c.amide('E394')


# #### chiral
# Side chains of Thr and Ile are chiral, incorrect atom labelling lead to the wrong chirality.  
# To fix use chiral('All | None | residue_list')

# In[22]:


st_c.chiral()


# #### Backbone
# Detects and fixes several problems with the backbone
# use any of 
# --fix_atoms All|None|Residue List 
# --fix_chain All|None|Break list
# --add_caps All|None|Terms|Breaks|Residue list
# --no_recheck
# --no_check_clashes
# 

# In[23]:


st_c.backbone()


# In[24]:


st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')


# #### fixside
# Detects and re-built missing protein side chains.   
# To fix use fixside('All | None | residue_list')

# In[25]:


st_c.fixside()


# #### getss
# Detects possible -S-S- bonds based on distance criteria.
# Proper simulation requires those bonds to be correctly set. Use All|None|residueList to mark them

# In[26]:


st_c.getss()


# In[27]:


st_c.getss('all')


# #### Add_hydrogens
#  Add Hydrogen Atoms. Auto: std changes at pH 7.0. His->Hie. pH: set pH value
#     list: Explicit list as [*:]HisXXHid, Interactive[_his]: Prompts for all selectable residues
#     Fixes missing side chain atoms unless --no_fix_side is set
#     Existing hydrogen atoms are removed before adding new ones unless --keep_h set.

# In[28]:


st_c.add_hydrogen()


# In[29]:


st_c.add_hydrogen('auto')


# #### clashes
# Detects steric clashes based on distance criteria.  
# Contacts are classified in: 
# * Severe: Too close atoms, usually indicating superimposed structures or badly modelled regions. Should be fixed.
# * Apolar: Vdw colissions.Usually fixed during the simulation.
# * Polar and ionic. Usually indicate wrong side chain conformations. Usually fixed during the simulation
# 

# In[30]:


st_c.clashes()


# Complete check in a single method

# In[31]:


st_c.checkall()


# In[32]:


st_c._save_structure(args['output_structure_path'])


# In[33]:


st_c.rem_hydrogen('yes')


# In[39]:


#st_c.add_hydrogen('--add_charges --add_mode auto')
#Alternative way calling through command line
import os
os.system('check_structure -i ' + args['output_structure_path'] + ' -o ' + args['output_structure_path_charges'] + ' add_hydrogen --add_charges --add_mode auto')


# In[35]:


#st_c._save_structure(args['output_structure_path_charges'])


# In[36]:


#st_c.revert_changes()


# In[ ]:





# In[ ]:




