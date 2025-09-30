# baktfold
Rapid &amp; standardized annotation of bacterial genomes, MAGs &amp; plasmids using protein structural information

## Install

```bash
conda create -n baktfold foldseek
conda activate baktfold
pip install -e .
baktfold --help
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 1 -d ~/baktfold_db/ 
```

* Where the `baktfold_db` for now is the Phold DB (for ProstT5)

```bash
foldseek databases Alphafold/Swiss-Prot swissprot tmp --threads 8
```
