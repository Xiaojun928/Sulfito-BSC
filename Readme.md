## Sulfitobacter speciation
Source code to replicate results for the manuscript "A Novel Bacterial Speciation Process Found in a Marine Symbiotic Population".

In the `01-phylogeny-construction` directory:
- Input:`.pep` and `.gene` files, the former contains amino acid sequences for each genome, and the latter contains nucleotide sequences for each genome;
- `S1.orthoFinder_blast.sh` and `S2.extract_scp_seq.pl` prepare the single copy core orthologs sequences in both amino acids (`.faa`) and nucleotides (`.dna`);
- Phylogenies based on  concatenated single copy orthologs alignments can be constructed using scripts `S3*-S4*`;
- Phylogenies of each single copy orthologs can be constructed using the scripts in  `nuc_seq`; 



In the `02-RMS-identification` directory:
- Scripts `s1*-s3*` identify the best hit restriction or modification enzymes for query genes;
- The reference (R-M system) data is downloaded form [REBASE]( http://rebase.neb.com/rebase/rebase.ftp.html);
- The `.pep` containing amino acid sequences for each genome is used as the query.

Note that the R-M system for each genome should be manually check whether the genes are neighbors after getting the best hits.


In the `03-Quartet-analyses` directory:
The number of supported topologies of each quatet was counted using [QuartetScores](https://github.com/lutteropp/QuartetScores)
The original QuatetScores did not output the count of each topology.
By modifying Ln117-120 in /your_installation_path/QuartetScores/src/QuartetCountConverter.hpp, and recompling the QuateteScores, the count of supported topologies for each quartet can be print.

Here, modified QuartetCountConverter.hpp was uploaded.

In the `04-identity-scanning` directory:
- Scripts `s1*-s2*` calculated the nucleotide identity between clade members and visualized results
- `s3.check_ANI.py` calculated the avg. ANI between clades


In the `05-neutral-simulation` directory:
