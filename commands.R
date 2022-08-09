library(pheatmap)
library(reticulate)
library(ggplot2)
library(dplyr)
library(jsonlite)
library(bio3d)

# load alphafold model

pd<-import("pandas")
features<-pd$read_pickle("full_outputs/alphafold/features.pkl")

str(features, max.level = 1)

models<-fromJSON("full_outputs/alphafold/ranking_debug.json")
models

best_model=pd$read_pickle("full_outputs/alphafold/result_model_2.pkl")
names(best_model)

best_model_plddt<-data.frame(residue=1:length(best_model$plddt), 
  plddt=best_model$plddt)

ggplot(best_model_plddt, aes(x=residue, y=plddt))+geom_point()

best_model_ptm<-pd$read_pickle("full_outputs/alphafold/result_model_3_ptm.pkl")
pheatmap(best_model_ptm$predicted_aligned_error, cluster_cols=F, cluster_rows=F)

# load rosetta model
np<-import("numpy")

rosetta_params<-np$load("full_outputs/rosettafold/t000_.3track.npz", mmap_mode="r")
rosetta_params$files
rosetta_params$f["dist"]

# multimer outputs

params=pd$read_pickle("multimer/term/result_model_1_multimer.pkl")
features=pd.read_pickle("multimer/term/features.pkl")

multi_plddt<-data.frame(residue=1:length(params$plddt), 
  plddt=params$plddt, chain=features$entity_id)

multi_plddt$protein<-ifelse(multi_plddt$chain==1, "eRF1", "eRF3")

ggplot(multi_plddt, aes(x=residue, y=plddt, color=protein))+geom_point()

pheatmap(params$predicted_aligned_error, cluster_cols=F, cluster_rows=F)

# load pdb plot beta factors
pdb<-read.pdb("multimer/term/ranked_0.pdb")

ggplot(unique(pdb$atom[, c("resno", "b", "chain")]), aes(x=resno, y=b, color=chain))+geom_point()+
  facet_grid(.~chain)


# calculate contacts

inds <- atom.select(pdb, "calpha")
contacts <- cmap(pdb, inds=atom.select(pdb, "calpha"), dcut=10, binary=F, collapse=F, mask.lower=F)

contacts_long<-reshape2::melt(contacts)

#select erf1 and erf3 residues
contacts_filt<-contacts_long %>%
  filter(Var1 %in% multi_plddt$residue[multi_plddt$protein=="eRF1"]) %>%
  filter(Var2 %in% multi_plddt$residue[multi_plddt$protein=="eRF3"]) %>%
  filter(value==1)


erf1_contacts<-unique(contacts_filt$Var1)
erf3_contacts<-unique(contacts_filt$Var2)

# separate into domains

pdb<-read.pdb("monomer/multi_domain/model_1.crderr.pdb")
domains<-read.table("monomer/multi_domain/DomainSpreadsheet-Homo_sapiens_Transcript_Domains_ENST00000545968.csv", 
  header=T, sep=",") %>%
  filter(Domain.source=="Gene3D")



for(i in 1:nrow(domains)){
  inds<-atom.selecT(pdb, resno=domains$Start[i]:domains$End[i])
  domain_atoms<-trim.pdb(pdb, inds)
  write.pdb(domain_atoms, file=paste0("domain_", i, ".pdb"))
}












