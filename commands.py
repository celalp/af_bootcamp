import pandas as pd
from Bio.PDB import *
import numpy as np
import json
import seaborn as sns
import pickle


# monomer
## read af output 

features=pd.read_pickle("full_outputs/alphafold/features.pkl")
features.keys()
features["template_sequence"]
features["template_domain_names"]

## model rankings
with open('full_outputs/alphafold/ranking_debug.json') as json_file:
    data = json.load(json_file)
    
data

best_model=pd.read_pickle("full_outputs/alphafold/result_model_2.pkl")

best_model.keys()
best_model["plddt"]
best_model_plddt=pd.DataFrame({"residue":range(len(best_model["plddt"])), "plddt":best_model["plddt"]})
sns.scatterplot(data=best_model_plddt, x="residue", y="plddt")


## load best model with ptm
best_ptm=pd.read_pickle("full_outputs/alphafold/result_model_3_ptm.pkl")
sns.heatmap(best_ptm["predicted_aligned_error"])


## rosettafold outputs
rosetta_params=np.load("full_outputs/rosettafold/t000_.3track.npz", mmap_mode="r")

for key in rosetta_params.keys():
    print(key)

# multimer outputs

with open("multimer/term/result_model_1_multimer.pkl", "rb") as pick:
    params=pickle.load(pick)

params.keys()

features=pd.read_pickle("multimer/term/features.pkl")

multi_plddt=pd.DataFrame({"residue":range(len(params["plddt"])), "plddt":params["plddt"], "chain":features["entity_id"]})
multi_plddt["protein"]="eRF1"
multi_plddt.loc[multi_plddt["chain"]==2, "protein"]="eRF3"
sns.scatterplot(data=multi_plddt, x="residue", y="plddt", hue="protein")

sns.heatmap(data=params["predicted_aligned_error"])

# calculate contacts

parser = PDBParser()
erf1_erf3=parser.get_structure(id="erf1_erf3", file="multimer/term/ranked_0.pdb")
chains=Selection.unfold_entities(erf1_erf3, "C")
chains

erf1=chains[0]
erf3=chains[1]

backbone_atoms=[]
for residue in erf1:
    backbone_atoms.append(residue["CA"])
    
ns=NeighborSearch(backbone_atoms)

erf3_contacts=[]

for residue in erf3:
    backbone_atom=residue["CA"]
    close_atoms=ns.search(backbone_atom.coord, 10)
    if len(close_atoms) > 0:
        erf3_residues=[atom.get_parent() for atom in close_atoms]
        erf3_contacts=erf3_contacts+erf3_residues
    

erf3_contacts=list(set(erf3_contacts))
erf3_contacts


erf1_abce1=parser.get_structure(id="erf1_erf3", file="multimer/recyc/ranked_0.pdb")


chains=Selection.unfold_entities(erf1_abce1, "C")

erf1=chains[0]
abce1=chains[1]

backbone_atoms=[]
for residue in erf1:
    backbone_atoms.append(residue["CA"])
    
ns=NeighborSearch(backbone_atoms)


abce1_contacts=[]

for residue in abce1:
    backbone_atom=residue["CA"]
    close_atoms=ns.search(backbone_atom.coord, 10)
    if len(close_atoms) > 0:
        abce1_residues=[atom.get_parent() for atom in close_atoms]
        abce1_contacts=abce1_contacts+abce1_residues
    

abce1_contacts=list(set(abce1_contacts))
abce1_contacts

same_residues=[residue for residue in abce1_contacts if residue in erf3_contacts]
same_residues.sort()
same_residues


# separate pdb into domans
parser = PDBParser()
wt=parser.get_structure(id="wt", file="monomer/multi_domain/model_1.crderr.pdb")

domains=pd.read_csv("monomer/multi_domain/DomainSpreadsheet-Homo_sapiens_Transcript_Domains_ENST00000545968.csv")
gene3d=domains[domains["Domain source"]=="Gene3D"]

starts=gene3d["Start"].tolist()
ends=gene3d["End"].tolist()

len(starts)==len(ends)
gene3d=domains[domains["Domain source"]=="Gene3D"]

class GetDomain(Select):
    def __init__(self, start, end):
        self.start=start
        self.end=end
        
    def accept_residue(self, residue):
        if residue.id[1] >= self.start and residue.id[1]<=self.end:
            return 1
        else:
            return 0

for i in range(len(starts)):
    filename="domain_"+str(i)+".pdb"
    io=PDBIO()
    io.set_structure(wt)
    io.save(filename, GetDomain(start=starts[i],end=ends[i]))









