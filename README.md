# baktfold
Rapid &amp; standardized annotation of bacterial genomes, MAGs &amp; plasmids using protein structural information

## Install

* We will be making baktfold available via pypi and conda shortly
* For now, please install from source
* The only non-python dependency is Foldseek

```bash
conda create -n baktfold foldseek
conda activate baktfold
pip install -e .
baktfold --help
```

* If you have a Mac with M-series Apple Silicon, you may need to install a particular version of Pytorch to utilise GPU-acceleration
* The same is true if you use other non-NVIDIA e.g. AMD GPUs
* See [this link](https://phold.readthedocs.io/en/latest/install/) for some more detail and further links

* Database installation

```bash
baktfold install --help
```

## Usage

* It is recommend you run baktfold with a GPU if you can.
* If you do not have a GPU, baktfold will still run, but the ProstT5 will be fairly slow.
* If you have a NVIDIA GPU, you can also use the `--foldseek-gpu` parameter to accelerate Foldseek

```bash
Usage: baktfold [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation          Print the citation(s) for this tool
  compare           Runs Foldseek vs baktfold db
  install           Installs ProstT5 model and baktfold database
  predict           Uses ProstT5 to predict 3Di tokens - GPU recommended
  proteins          baktfold protein-predict then comapare all in one - GPU...
  proteins-compare  Runs Foldseek vs baktfold db on proteins input
  proteins-predict  Runs ProstT5 on a multiFASTA input - GPU recommended
  run               baktfold predict then comapare all in one - GPU...
```

* The two most useful commands are `baktfold run` and `baktfold proteins`
* `baktfold run` accepts a Bakta json file as its input, and by default, it will annotate all hypothetical CDS and return a variety of Bakta-like compliant output formats
* `baktfold proteins` accepts a protein FASTA `.faa` format file as input. It will annotate all protein sequences and return a variety of `bakta_proteins`-like output formats
* `baktfold predict` and `baktfold compare` split `baktfold run` into the ProstT5 and Foldseek modules, while `baktfold proteins-predict` and `baktfold proteins-compare` do the same for `baktfold proteins` (useful if you have non-NVIDIA GPUs)


## Example

```bash
# with nvidia gpu 
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   --foldseek-gpu
# without nvidia gpu available
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   
```
