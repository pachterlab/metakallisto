These python scripts allow metagenomic sequence data to be analyzed with the fast, accurate RNA-Seq abundance estimator [kallisto](https://github.com/pachterlab/kallisto). Both taxa identification and abundance estimation can be performed at the exact-genome level, as demonstrated in our paper ["Pseudoalignment for metagenomic read assignment"](https://arxiv.org/abs/1510.07371).

##Prerequisites

The full pipeline as shown in our paper requires the installation of genome distance estimator [Mash](https://github.com/marbl/mash), RNA-Seq abundance estimator [kallisto](https://github.com/pachterlab/kallisto), and python 2.7+. 

Python modules needed for running the scripts include [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), and [Biopython](http://biopython.org/wiki/Download). 

##Usage

`compare_metagenomic_results_to_truth.py` compares the output of a range of metagenomic analysis tools such as kallisto to the ground truth of the [illumina 100](http://dx.doi.org/10.1371/journal.pone.0031386) metagenomic dataset used in our paper. It is run from the command line:
```
python compare_metagenomic_results_to_truth.py [filename] [program (default kallisto)] 
						<-g show graphs?> <-s save graphs?> <--taxa level (defaults to all)>

positional arguments:
  filename              Output file of metagenomic analysis tool
  program               Source program that created output. Valid options are:
                        kallisto, kraken, clark, gasic, express. Defaults to
                        kallisto.

optional arguments:
  -h, --help            show this help message and exit
  -g, --show-graphs     Display graphs of calculated errors
  -s, --save-graphs     Save graphs of calculated errors to file
  --taxa TAXA           Desired taxa level of analysis. Accepts one of:
                        strain, species, genus, phylum. Defaults to all
                        levels.
  --dataset DATASET     Dataset truth to be compared to. Accepts: i100,
                        no_truth. Defaults to i100.
  --bootstraps BOOTSTRAPS
                        Directory containing .tsv files for kallisto
                        bootstraps, to be converted into errors.
```
`plotfunctions.py` is a helper module necessary to allow compare_metagenomic_results_to_truth.py to plot the results of analysis.

`mash_kallisto_pipeline.py` will process the output of Mash (run on a large set of metagenomic reference genomes), and select the top N (user-defined) genomes that match each species that Mash identified in the raw sequenced reads. Python script should be run from the directory containing the raw .fa or .mfa files, and will move matching files to a user-specified directory for later indexing. The script will also process each genome to concatenate contigs/chromosomes, and attempt to look up taxids based on common UIDs present in the file name; a taxid will allow metagenomic programs such as Kraken or CLARK to use these genomes to make a custom database.

```
mash_kallisto_pipeline.py [-h] [--directory DIRECTORY] [--dry-run] filename top_strains

positional arguments:
  filename              Mash output file
  top_strains           How many strains of each species to keep for the
                        quantification step

optional arguments:
  -h, --help            show this help message and exit
  --directory DIRECTORY
                        Directory to put files for kallisto index creation.
                        Default is moving them one directory up.
  --dry-run             If set, lists files that would be moved, but does not
                        create or move any files.
```
