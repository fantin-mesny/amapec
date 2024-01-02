from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import math
from math import sqrt, exp
import argparse
import sys
import os
import numpy as np
import subprocess
import multiprocessing
#import warnings
#warnings.filterwarnings(action='ignore', category=PDBConstructionWarning)#, module='sklearn')

def get_params(argv):
	parser = argparse.ArgumentParser(description='Calculate biochemical properties from structures')
	parser.add_argument('-i', '--i', help="Comma-separated directories including structures'", required=True)
	parser.add_argument('-o', '--o', help="Output/working directory'", required=True)
	parser.add_argument('-nproc', '--nproc', help="Number of threads", default=100)
	a = parser.parse_args()
	return a
	

def get_mean_pLDDT(File): #optimized on AF2 structures
	with open(File,'r') as inp:
		lines=inp.readlines()
	lines_list=[l.split(' ') for l in lines]
	lines_list_clean=[[e for e in l if e!=''] for l in lines_list if l[0]=='ATOM']
	lines_list_filtered=[(str(l[4].replace('A','A,')+','+l[5]).replace(',,',',').split(',')[1],l[-3]) for l in lines_list_clean]
	vals=pd.DataFrame(lines_list_filtered)
	vals=vals.drop_duplicates()
	return vals[1].astype(float).mean()

######################################################################################################
############### Run and parse Fpocket ##################################
# user manual: https://fpocket.sourceforge.net/manual_fpocket2.pdf

def runAndParseFpocket(File):
	File=File.replace('|','\|')
	outp=subprocess.check_output(['fpocket', '-f',File,'-d'])
	outp2=str(outp)
	subprocess.check_output(['rm','-r',File.replace('.pdb','_out')])
	outp2=outp2.split('\\n')
	return pd.DataFrame([line.split(' ')[:-3] for line in outp2][1:-1],columns=outp2[0].split(' ')[:-3],dtype=float)

def getPocketInfo(File):
	pockets=runAndParseFpocket(File)
	pocketsDF=pd.DataFrame(columns=['value'])
	pocketsDF.loc['Number of pockets predicted','value']=len(pockets)
	pocketsDF.loc['Mean pocket volume','value']=pockets['volume'].mean()
	pocketsDF.loc['Volume of the biggest pocket','value']=pockets['volume'].max()
	pocketsDF.loc['Volume of the smallest pocket','value']=pockets['volume'].min()
	pocketsDF.loc['Mean number of alpha spheres per pocket','value']=pockets['nb_asph'].mean()
	pocketsDF.loc['Mean density of alpha spheres per pocket','value']=pockets['as_density'].mean()
	pocketsDF.loc['Mean flexibility of pockets','value']=pockets['flex'].mean()
	pocketsDF.loc['Mean hydrophobicity score per pocket','value']=pockets['hydrophobicity_score'].mean()
	pocketsDF.loc['Mean polarity score per pocket','value']=pockets['polarity_score'].mean()
	pocketsDF.loc['Mean charge score per pocket','value']=pockets['charge_score'].mean()
	pocketsDF.loc['Mean volume score per pocket','value']=pockets['volume_score'].mean()
	return pocketsDF.fillna(0)


######################################################################################################
############### Taken from https://www.kaggle.com/code/saidtaghadouini/exploring-the-3d-geometry#Sphericity
def compute_sphericity(File):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', File)
    
    # Get the coordinates of the atoms in the structure
    coords = []
    for residue in structure.get_residues():
        # Calculate the center of the amino acid as the mean of the coordinates of its atoms
        center = np.mean([atom.get_coord() for atom in residue], axis=0)
        coords.append(center)
    
    # Convert the coordinates to a NumPy array
    coords = np.array(coords)
    
    # Compute the sphericity
    I = np.zeros(3)
    for i in range(3):
        I[i] = (np.sum((coords[:,i]-coords[:,i].mean())**2))/len(coords)
    sphericity = (3*np.sqrt(I[0]*I[1]+I[1]*I[2]+I[2]*I[0])/(2*(I[0]+I[1]+I[2])))
    
    return sphericity




######################################################################################################
################ Taken from ssbio https://github.com/SBRG/ssbio/blob/08c1ce375a8ee6a4b77c5d2fdc93fcb5e5fe97c0/ssbio/protein/structure/properties/residues.py#L24
	
	
def search_ss_bonds(model, threshold=3.0):
    """ Searches S-S bonds based on distances
        between atoms in the structure (first model only).
        Average distance is 2.05A. Threshold is 3A default.
        Returns iterator with tuples of residues.
        ADAPTED FROM JOAO RODRIGUES' BIOPYTHON GSOC PROJECT (http://biopython.org/wiki/GSOC2010_Joao)
    """

    # Taken from http://docs.python.org/library/itertools.html
    # Python 2.4 does not include itertools.combinations

    def combinations(iterable, r):
        # combinations('ABCD', 2) --> AB AC AD BC BD CD
        # combinations(range(4), 3) --> 012 013 023 123
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = list(range(r))
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i + 1, r):
                indices[j] = indices[j - 1] + 1
            yield tuple(pool[i] for i in indices)

    cysteines = [r for r in model.get_residues() if r.get_resname() == 'CYS']

    pairs = combinations(cysteines, 2)  # Iterator with pairs

    bridges = []
    for cys_pair in pairs:
        try:
            if cys_pair[0]['SG'] - cys_pair[1]['SG'] < threshold:
                bridges.append(cys_pair)
        except KeyError:  # This will occur when a CYS residue is missing a SG atom for some reason
            log.error('{}: no SG atom found for one or both of the cysteine residues {}'.format(model, cys_pair))
            continue

    return bridges
	

########################################################################################
#### This part is stolen from https://github.com/ugSUBMARINE/structural-properties #####

sasa = {"ALA": 75,"CYS": 115,"ASP": 130,"GLU": 161,"PHE": 209,"GLY": 0,"HIS": 180,"ILE": 172,"LYS": 205,"LEU": 172,"MET": 184,"ASN": 142,"PRO": 134,"GLN": 173,"ARG": 236,"SER": 95,"THR": 130,"VAL": 143,"TRP": 254,"TYR": 222} # surface accessible side chain area

def data_coord_extraction(
    target_pdb_file: str,
) -> tuple[
    np.ndarray[tuple[int, int], np.dtype[str]],
    np.ndarray[tuple[int, int], np.dtype[float]],
]:
    """extracts the coordinates and the residue data from a pdb file
    :parameter
         - target_pdb_file:
           path to pdb file for protein of interest
    :returns
         - new_data: 2D ndarray
           contains information about all residues like [[Atom type,
           Residue 3letter, ChainID, ResidueID],...]
         - new_coords:
           contains coordinates of corresponding residues to the new_data
           entries
    """
    # list of all data of the entries like [[Atom type, Residue 3letter,
    # ChainID, ResidueID],...]
    res_data = []
    # list of all coordinates of the entries like [[x1, y1, z1],...]
    res_coords = []
    # reading the pdb file
    file = open(target_pdb_file, "r")
    for line in file:
        if "ATOM  " in line[:6]:
            line = line.strip()
            res_data.append(
                [
                    line[12:16].replace(" ", "").strip(),
                    line[17:20].replace(" ", "").strip(),
                    line[21].replace(" ", "").strip(),
                    line[22:26].replace(" ", "").strip(),
                ]
            )
            res_coords.append(
                [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]
            )
    file.close()

    res_data = np.asarray(res_data)
    res_coords = np.asarray(res_coords, dtype=float)
    return res_data, res_coords


def dist_calc(
    arr1: np.ndarray[tuple[int, int], np.dtype[int | float]],
    arr2: np.ndarray[tuple[int, int], np.dtype[int | float]],
) -> np.ndarray[tuple[int, int], np.dtype[float]]:
    """calculates distance between arr1 and arr2 and returns a 2D array with
    all distances of all arr1 points against all arr2 points
     :parameter
         - arr1, arr2:
           2D arrays of 1D lists with 3D coordinates eg [[x1, y1, z1],...]
     :return
         - dist:
           len(arr1) x len(arr2) distance matrix between arr1 and arr2"""
    # get only the x,y,z coordinates from the input arrays and reshape them,
    # so they can be subtracted from each other
    arr1_coords_rs = arr1.reshape(arr1.shape[0], 1, arr1.shape[1])
    arr2_coord_rs = arr2.reshape(1, arr2.shape[0], arr1.shape[1])
    # calculating the distance between each point and returning a 2D array
    # with all distances
    dist = np.sqrt(((arr1_coords_rs - arr2_coord_rs) ** 2).sum(axis=2))
    return dist


def c_cluster(pair_num: list[list[str]]) -> list[list[str]]:
    """clusters pairs of interactions that feature on common residue
    :parameter
        - pair_num:
          pairs of salt bridge forming residues as their residue number
    :return
        - pair_num:
          all residues that are in the same cluster in one list

    """
    pair_num = pair_num.tolist()
    # which residues are forming salt bridges and how often they are in contact
    # with another residue
    multi_ap, mac = np.unique(pair_num, return_counts=True)
    # dict describing the number of contacts for each residue
    app_dict = dict(zip(multi_ap, mac))

    """
    while the number of clusters changes due to merging
       for all existing clusters (starts with pairs)
            if the residues in the current cluster appear only in one cluster
                put this cluster in the clusters list
            else
                check for each cluster if one of the residues of the current
                cluster that's looked at is in another cluster
                if that's the case
                    put the current cluster and the one where one of the
                    residues are in in the inter_cluster
                make the residues in inter_cluster unique
                if this cluster is not in the clusters list put it in there
            change the clusters that are looked at to the new created clusters
            list and restart if the number of clusters is not the same as in
            the previous round
    """
    prev_len = 0
    while True:
        clusters = []
        for i in pair_num:
            if np.sum(list(map(app_dict.get, i))) == len(i):
                clusters.append(i)
            else:
                inter_cluster = []
                for k in pair_num:
                    for j in i:
                        if j in k:
                            inter_cluster += i
                            inter_cluster += k
                potential_cluster = np.unique(inter_cluster).tolist()
                if potential_cluster not in clusters:
                    clusters.append(potential_cluster)
        pair_num = clusters
        c_len = len(clusters)
        if c_len == prev_len:
            break
        else:
            prev_len = c_len
    return pair_num


def join_res_data(
    part0: np.ndarray[tuple[int, int], np.dtype[str]],
    part1: np.ndarray[tuple[int, int], np.dtype[str]],
) -> np.ndarray[tuple[int, int], np.dtype[str]]:
    joined_data = []
    for i in range(len(part0)):
        joined_data.append(["-".join(part0[i]), "-".join(part1[i])])
    return np.asarray(joined_data)


def create_output(
    pair_num: list[list[str]],
    pair_num_ori: list[list[str]],
    create_file: str = None,
    silent: bool = False,
) -> None:
    """prints output in the terminal and optionally creates a output file
    :parameter
        - pair_num:
          clustered residues in one list
        - pair_num_ori:
          the original pairs of interacting residues
        - create_file:
          how the output file will be named
        - silent:
          whether to print output in terminal or not
    :return
        - None
    """
    # list with cluster sizes
    cluster_sizes = [len(i) for i in pair_num]
    # contacts per cluster_sizes
    cpc = []
    for i in pair_num:
        inter_pairs = []
        for k in i:
            inter_pairs += pair_num_ori[np.where(pair_num_ori == k)[0]].tolist()
        cpc.append(np.unique(inter_pairs, axis=0).shape[0])
        r = []
        for f in np.unique(inter_pairs):
            f_split = f.split("-")
            r.append(f_split[-1])
    # output
    output=[]
    if create_file is not None:
        data_file = open(f"{create_file}.csv", "w+")
        data_file.write("InteractingResidues,ContactsPerCluster\n")
    for i in range(len(pair_num)):
        output.append({'Interacting Residues':" - ".join(pair_num[i]), 'Contacts Per Cluster':cpc[i]})
    return pd.DataFrame(output)


def find_saltbridges(
    file_path: str,
    max_dist: int | float = 3.5,
    sele_chain: str = None,
    create_file: str = None,
    silent: bool = False,
):
    """searches for saltbridges in a given structure with the option to limit the chains
    :parameter
        - file_path:
          path to the pdb file
        - sele_chain:
          to select a specific chain use e.g. 'A' in which the salt bridges should be
          calculated
        - create_file:
          how the output file will be named - gets split if '_' present
        - silent:
          whether to print output in terminal or not
    :return
        - None
    """
    data, coords = data_coord_extraction(file_path)
    # ARG NH1, NH2
    # LYS NZ
    # HIS ND1 NE2
    # ASP OD1 OD2
    # GLU OE1 OE2
    aa = ["ARG", "LYS", "HIS", "ASP", "GLU"]
    atom = ["NH1", "NH2", "NZ", "ND1", "NE2", "OD1", "OD2", "OE1", "OE2"]
    charge = {"ARG": 1, "LYS": 1, "HIS": 1, "ASP": -1, "GLU": -1}
    tests = []
    for i in aa:
        tests.append((data[:, 1] == i).tolist())
    for i in atom:
        tests.append((data[:, 0] == i).tolist())
    tests = np.asarray(tests)
    # which data entry contains the right amino acid and the right atom type
    test_conf = np.sum(tests, axis=0) == 2
    if sele_chain is None:
        chain_test = np.ones(data[test_conf].shape[0]).astype(bool)
    else:
        # to get data entries from selected chain(s)
        chain_test = data[test_conf][:, 2] == sele_chain
    # ResidueIDs for interacting residues that are able to form bridges
    sele_data = data[test_conf][chain_test]
    # their coordinates
    sele_coords = coords[test_conf][chain_test]
    # distance matrix between all the potential residues
    dists = dist_calc(sele_coords, sele_coords)
    # all atoms in salt bridge distance
    dists = dists < max_dist
    # to only get interactions once
    dists = np.triu(dists, 1)
    # indices of which sele_data are in interacting distance
    pairs = np.where(dists)
    # charge of the interacting residues
    charge_map = np.asarray(list(map(charge.get, sele_data[:, 1])))
    # only allow interacting pairs to be valid if their charge is opposite
    valid_charge_bool = charge_map[pairs[0]] + charge_map[pairs[1]] == 0
    # ResidueIDs of the valid interaction pairs
    valid_pairs = np.column_stack(
        (sele_data[:, [1, 2, 3]][pairs[0]], sele_data[:, [1, 2, 3]][pairs[1]])
    )[valid_charge_bool]
    valid_pairs = np.unique(valid_pairs, axis=0)

    # pairs of salt bridge forming residues as their residue number as strings
    pair_num = join_res_data(valid_pairs[:, 3:], valid_pairs[:, :3])
    pair_num_ori = pair_num

    pair_num = c_cluster(pair_num)

    return create_output(pair_num, pair_num_ori, create_file=create_file, silent=silent)
    #return pair_num_ori

def create_output_hy(
    pair_num: list[list[str]],
    pair_num_ori: list[list[str]],
    create_file: str = None,
    silent: bool = False,
) -> None:
    """prints output in the terminal and optionally creates a output file
    :parameter
        - pair_num:
          clustered residues in one list
        - pair_num_ori:
          the original pairs of interacting residues
        - create_file:
          how the output file will be named
        - silent:
          whether to print output in terminal or not
    :return
        - None
    """
    # list with cluster sizes
    cluster_sizes = [len(i) for i in pair_num]
    # contacts per cluster_sizes
    cpc = []
    # surface area per cluster
    spc = []
    for i in pair_num:
        inter_pairs = []
        for k in i:
            inter_pairs += pair_num_ori[np.where(pair_num_ori == k)[0]].tolist()
        cpc.append(np.unique(inter_pairs, axis=0).shape[0])
        r = []
        surface_area = 0
        for f in np.unique(inter_pairs):
            f_split = f.split("-")
            surface_area += int(sasa[f_split[0]])
            r.append(f_split[-1])
        spc.append(surface_area)
    output_hy=[]
    for i in range(len(pair_num)):
        output_hy.append({'Interacting Residues':" - ".join(pair_num[i]),'Contacts Per Cluster':cpc[i],'Surface Area Per Cluster':str(spc[i])})
    return pd.DataFrame(output_hy)


def hydr_cluster(
    file_path: str,
    sele_chain: str = None,
    create_file: str = None,
    silent: bool = False,
) -> None:
    """calculates hydrophobic cluster with the option to select a chain
    :parameter
        - file_path:
          path to the pdb file
        - sele_chain:
          to select a specific chain use e.g. 'A' in which the salt bridges
          should be calculated
        - create_file:
          how the output file will be named - gets split if '_' present
        - silent:
          whether to print output in terminal or not
    :return
        - None
    """
    data, coords = data_coord_extraction(file_path)
    aa = ["ILE", "LEU", "VAL"]
    atom = ["N", "H", "CA", "HA", "C", "O"]
    tests = []
    for i in aa:
        tests.append((data[:, 1] == i).tolist())
    for i in atom:
        tests.append((data[:, 0] != i).tolist())
    heavy_atom = []
    for i in data[:, 0]:
        heavy_atom.append(not i.startswith("H"))
    tests.append(heavy_atom)
    tests = np.asarray(tests)
    # which data entry contains the right amino acid and the right atom type
    test_conf = np.sum(tests, axis=0) == 8
    if sele_chain is None:
        chain_test = np.ones(data[test_conf].shape[0]).astype(bool)
    else:
        # to get data entries from selected chain(s)
        chain_test = data[test_conf][:, 2] == sele_chain
    # ResidueIDs for interacting residues that are able to form bridges
    sele_data = data[test_conf][chain_test]
    # their coordinates
    sele_coords = coords[test_conf][chain_test]
    # distance matrix between all the potential residues
    dists = dist_calc(sele_coords, sele_coords)
    # all atoms in hydrophobic interaction distance
    dists = dists < 6.56
    # to only get interactions once
    dists = np.triu(dists, 1)
    pair_ind0, pair_ind1 = np.where(dists)
    excl_same = np.any(
        sele_data[:, [1, 2, 3]][pair_ind0] != sele_data[:, [1, 2, 3]][pair_ind1], axis=1
    )
    # pairs of hydrophobic interactions with Amino Acid, Chain and Number
    valid_pairs = np.column_stack(
        (
            sele_data[:, [1, 2, 3]][pair_ind0][excl_same],
            sele_data[:, [1, 2, 3]][pair_ind1][excl_same],
        )
    )
    # to only have one entry per residues interaction and not of all their
    # atoms
    valid_pairs = np.unique(valid_pairs, axis=0)

    # pairs of cluster forming residues as their residue number as strings
    pair_num = join_res_data(valid_pairs[:, 3:], valid_pairs[:, :3])
    pair_num_ori = pair_num

    pair_num = c_cluster(pair_num)

    return create_output_hy(pair_num, pair_num_ori, create_file=create_file, silent=silent)




########################################################################################


def Rg(filename): # from https://github.com/sarisabban/Rg
	'''
	Calculates the Radius of Gyration (Rg) of a protein given its .pdb 
	structure file. Returns the Rg integer value in Angstrom.
	'''
	coord = list()
	mass = list()
	Structure = open(filename, 'r')
	for line in Structure:
		try:
			line = line.split()
			x = float(line[6])
			y = float(line[7])
			z = float(line[8])
			coord.append([x, y, z])
			if line[-1] == 'C':
				mass.append(12.0107)
			elif line[-1] == 'O':
				mass.append(15.9994)
			elif line[-1] == 'N':
				mass.append(14.0067)
			elif line[-1] == 'S':
				mass.append(32.065)
		except:
			pass
	xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
	tmass = sum(mass)
	rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
	mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
	rg = math.sqrt(rr / tmass-mm)
	return(round(rg, 3))

def readPDB(pdb_file):
    """
    Takes a pdb file and reads in the coordinates of each titratable group.
    Assigns pka and charge state of each group.
    """

    TITR_ATOM = {"ASP":"CG ",
                 "GLU":"CD ",
                 "TYR":"OH ",
                 "ARG":"CZ ",
                 "HIS":"NE2",
                 "LYS":"NZ "}

    PKA_DICT = {"ASP":     4.0,
                "CYS":     8.5,
                "GLU":     4.4,
                "TYR":    10.0,
                "CTERM":   3.1,
                "ARG":    12.0,
                "HIS":     6.5,
                "LYS":    10.4,
                "NTERM":   8.0}

    CHARGE_DICT = {"ASP":   -1.,
                   "CYS":   -1.,
                   "GLU":   -1.,
                   "TYR":   -1.,
                   "CTERM": -1.,
                   "ARG":    1.,
                   "HIS":    1.,
                   "LYS":    1.,
                   "NTERM":  1.}

    # Open pdb_file and read each line into pdb (a list of lines)
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Grab only ATOM entries that are titratable
    pdb = [l for l in pdb if l[0:4] == "ATOM" and
                             l[17:20] in list(TITR_ATOM.keys()) and
                             l[13:16] == TITR_ATOM[l[17:20]]]

    # Initalize lists to hold coordinates, pkas, and charge
    coord, pKa, charge = [], [], []

    # Go through each line in the pdb file
    for line in pdb:
        amino_acid = line[17:20]

        # Grab the xyz coordinates
        coord.append([float(line[30+8*i:38+8*i]) for i in range(3)])

        pKa.append(PKA_DICT[amino_acid])

        # Look up charge
        charge.append(CHARGE_DICT[amino_acid])

    # Return the coordinates, pka, and charge
    return coord, pKa, charge

def hendersonHasselbach(pKa,charge,pH):
    """
    Calculate the fractional charge on a group with pKa and charge at some
    pH value.
    """

    return charge/(1 + 10**(charge*(pH-pKa)))


def pdbCoulomb(coord,pKa,charge,dielec_const,ionic_str,pH,temperature):
    """
    Calculates the energy of a structure given the coordinates of each
    charged atom, their fractional charge, the dielectric constant, and the
    ionic strength.
    """

    ionic_str = ionic_str/1000

    # Initialize variables
    kappa = 50.29*sqrt(ionic_str/(dielec_const*temperature))
    num_groups = len(coord)
    energy = 0.

    hh_chg = [hendersonHasselbach(pKa[i],charge[i],pH)
              for i in range(num_groups)]

    # Calculate energy of interaction of every ij interaction (making sure not
    # to double count; note we start j at i + 1).
    for i in range(num_groups):
        for j in range(i+1,num_groups):

            # Calculate distance between atom i and atom j
            r = sqrt(sum([(coord[i][k]-coord[j][k])**2 for k in range(3)]))

            # Add the energy of this interaction to the total energy
            energy += 332*hh_chg[i]*hh_chg[j]/(r*dielec_const)*exp(-kappa*r)

    # Return energy
    return energy
    
    
    
    
def getStructureProperties(Dir, struc,dfs):
			parser=PDBParser()
			structure=parser.get_structure('TEST', Dir+struc)
			model=structure[0]
			
			dssp=DSSP(model, Dir+struc, dssp='mkdssp')
			ss={'H':'Structure percentage in helices','G':'Structure percentage in helices','B':'Structure percentage in beta structures','E':'Structure percentage in beta structures','I':'Structure percentage in coils','T':'Structure percentage in coils','S':'Structure percentage in coils','-':'Structure percentage with undefined secondary structure'} # source: https://doi.org/10.1371/journal.pone.0028464
			renamed_dssp={'H':'Structure percentage in 4-turn helices','G':'Structure percentage in 3-turn helices','I':'Structure percentage in 5-turn helices', 'T':'Structure percentage in hydrogen bonded turns','E':'Structure percentage in beta sheets','B':'Structure percentage in residues in isolated beta bridges','S':'Structure percentage in bends','-':'-'}
			##"Along with the DSSP program, the helical-structure includes residues from α- and 310–helices. The beta-structure includes residues from isolated β-bridges and extended strands involved in β-sheets. Residues from π-helices, hydrogen-bonded turns and bends are included in the irregular structure (coil)."
			dssp=pd.DataFrame(dssp)
			dssp['SS']=dssp[2].map(ss)
			dssp[2]=dssp[2].map(renamed_dssp)
			df1=pd.DataFrame(dssp[2].value_counts()/len(dssp)).rename(columns={2:'value'})
			df2=pd.DataFrame(dssp['SS'].value_counts()/len(dssp)).rename(columns={'SS':'value'})
			df=pd.concat([df1,df2])
			for var in list(ss.values())+list(renamed_dssp.values()):
				if var not in df.index:
					df.loc[var,'value']=0
			df=df.drop(index=['-'])
			df.loc['Accessible surface area','value']=dssp[3].sum()
			df.loc['Hydrophobic accessible surface area','value']=dssp[dssp[1].isin(['A','V','L','I','M','F','P'])][3].sum()
			hydrophobicity={'G':0,'Q':0,'S':0.07,'T':0.07,'N':0.09,'D':0.66,'E':0.67,'R':0.85,'A':0.87,'H':0.87,'C':1.52,'K':1.64,'M':1.67,'V':1.87,'L':2.17,'Y':2.76,'P':2.77,'F':2.87,'I':3.15,'W':3.77} # values from https://doi.org/10.1186/1471-2105-7-S4-S14
			df.loc['Surface hydrophobicity']=sum([dssp.loc[ind,3]*hydrophobicity[dssp.loc[ind,1]] for ind in dssp.index])/dssp[3].sum()
			df.loc['Acidic accessible surface area','value']=dssp[dssp[1].isin(['D','E'])][3].sum()
			df.loc['Basic accessible surface area','value']=dssp[dssp[1].isin(['K','R','H'])][3].sum()
			df.loc['Polar accessible surface area','value']=dssp[dssp[1].isin(['N','Q','Y','S','T'])][3].sum()
			for aa in hydrophobicity:
				df.loc['Accessible surface area occupied by amino acid '+aa]=dssp[dssp[1]==aa][3].sum()
				
				
			df.loc['Radius of gyration','value']=Rg(Dir+struc)
			df.loc['Number of disulphide bonds','value']=len(search_ss_bonds(model))
			
			saltBridges=find_saltbridges(Dir+struc)
			if len(saltBridges)>0:
				df.loc['Number of salt bridges','value']=saltBridges['Contacts Per Cluster'].sum()
			else:
				df.loc['Number of salt bridges','value']=0
				
				
			df.loc['Sphericity','value']=compute_sphericity(Dir+struc)
				
			hydrophobicClusters=hydr_cluster(Dir+struc)
			if len(hydrophobicClusters)>0:
				df.loc['Surface area covered by hydrophobic clusters','value']=float(hydrophobicClusters['Surface Area Per Cluster'].sum())
			else:
				df.loc['Surface area covered by hydrophobic clusters','value']=0
			df.loc['Number of hydrophobic clusters','value']=len(hydrophobicClusters)

			# coulomb energy
			dielec_const=40.0
			ionic_str=100
			pH=7
			temperature=298
			coord, pKa, charge=readPDB(Dir+struc)
			coulomb=pdbCoulomb(coord,pKa,charge,dielec_const,ionic_str,pH,temperature)
			df.loc['Coulomb energy','value']=coulomb
			df=pd.concat([df,getPocketInfo(Dir+struc)]).rename(columns={'value':struc.replace('.pdb','')})
			dfs.append(df)
			#print(df)


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	pool = multiprocessing.Pool(int(a.nproc))
	manager = multiprocessing.Manager()
	dfs=manager.list()
	
	if 'ignored.txt' in os.listdir(a.o):
		with open(a.o+'/ignored.txt') as inp:
			inpLines=inp.readlines()
		toIgnore=[ID.replace('\n','.pdb') for ID in inpLines]
	else:
		toIgnore=[]


	for Dir in a.i.split(','):
		if not a.i.endswith('/'): 
			Dir=Dir+'/'
		for struc in [st for st in os.listdir(Dir) if st.endswith('.pdb')]:
			if struc not in toIgnore:
				pool.apply_async(getStructureProperties, args=(Dir,struc,dfs))
	
	pool.close()
	pool.join()

	DF=pd.concat(dfs,axis=1).T.fillna(0)

	DF.to_csv(a.o+'/strucProperties.csv')
