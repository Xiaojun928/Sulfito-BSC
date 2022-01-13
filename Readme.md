## Sulfitobacter speciation
Source code to replicate results for the manuscript "A Novel Bacterial Speciation Process Found in a Marine Symbiotic Population".

In the `01-phylogeny-construction` directory:
- Input:`.pep` and `.gene` files, the former contains amino acid sequences for each genome, and the latter contains nucleotide sequences for each genome;
- `S1.orthoFinder_blast.sh` and `S2.extract_scp_seq.pl` prepare the single copy core orthologs sequences in both amino acids (`.faa`) and nucleotides (`.dna`);
- Phylogenies based on 16S rRNA genes or concatenated single copy orthologs alignments can be constructed using scripts `16S_tree_constr.sh` and `S3*-S4*`;
- Phylogenies of each single copy orthologs can be constructed using the scripts in `aa_seq` and `nuc_seq`; 
- `topo_check_outgroup_rooted.R` is used to determined the topology supported by each gene tree in `./with_outgroup`(gene trees each was constructed using 28 members of three clades and two outgroups based on either amino acids or nucleotides).
- Gene trees in `./without_outgroup` each was constructed using only 28 members of three clades, and rooted by [MAD](http://www.nature.com/articles/s41559-017-0193) and [MV](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182238) methods, respectively.  The supporting topology of each MAD-rooted gene tree was determined using `topo_check_MADroot.R`, and MV-rooted gene tree was determined with `topo_check_MVroot.R`.


In the `02-RMS-identification` directory:
- Scripts `s1*-s3*` identify the best hit restriction or modification enzymes for query genes;
- The reference (R-M system) data is downloaded form [REBASE]( http://rebase.neb.com/rebase/rebase.ftp.html);
- The `.pep` containing amino acid sequences for each genome is used as the query.

Note that the R-M system for each genome should be manually check whether the genes are neighbors after getting the best hits.

