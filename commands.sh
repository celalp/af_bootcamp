bash

#lddt command
lddt -c -p stereo_chemical_props.txt full_outputs/alphafold/ranked_0.pdb full_outputs/rosettafold/model/model_1.crderr.pdb

# mustang command
cd /MUSTANG_v3.2.3/
make

cd /home
cd monomer/allergens
files=$(ls)
mustang-3.2.3 -i $files -r ON