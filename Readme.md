## Sulfitobacter speciation
Source code to replicate results for the manuscript "Speciation of an Episymbiotic Marine Bacterial Population Parallels that of Eukaryotes Under Incomplete Lineage Sorting".

In the `01-phylogeny-construction` directory:
- Input:`.faa` and `.gene` files, the former contains amino acid sequences for each genome, and the latter contains nucleotide sequences for each genome;
- `S1.orthoFinder_blast.sh` and `S2.extract_scp_seq.pl` prepare the singel copy core orthologs sequences in both amino acids and nucleotides;
- Phylogenies based on 16S rRNA genes or concatenated single copy orthologs alignments can be constructed using scripts `16S_tree_constr.sh` and `S3*-S4*`;
- Phylogenies of each single copy orthologs can be constructed using the scripts in `aa_seq` and `nuc_seq`; 
- `topo_determine.R` is used to determined the topology supported by each gene tree (based on either amino acids or nucleotides).

In the `02-demography-inference` directory:
- XXX is used to infer the demography of Sulfitobacter population
- ...(**TDB**)

In the `03-RMS-identification` directory:
- Scirpts `s1*-s3*` identify the best hit restriction or modification enzymes for query genes;
- The reference (R-M system) data is downloaded form [REBASE]( http://rebase.neb.com/rebase/rebase.ftp.html);
- The `.faa` containing amino acid sequences for each genome is used as the query.

Note that the R-M system for each genome should be manually check whether the genes are neighbors after getting the best hits.


