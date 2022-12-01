# Nanopore xGen Amplicon Assembly Pipeline 
This pipeline automates the process of assembling genomes from IDT's xGen Amplicon Panels sequenced using Oxford Nanopore technology.

## Technical Considerations

### xGen Primer File

For primer clipping, this pipeline uses the tool [Primerclip](https://github.com/swiftbiosciences/primerclip) which was developed for IDT's xGen technology (formerly Swift Biosciences) specifically. 

The input for Primerclip is **different** from the traditional BED file used for tools such as iVar. Instead, IDT creates a special "Masterfile" text file that contains information about each potential amplicon (amplicon coordinates, primer coordinates, primer sequences, etc). IDT will provided this file to those who have purchased one of their products.

This file is requires for the pipeline, and must be specified via the ```--primers``` option.
### Medaka Model

Medaka requires information about the pore, sequencing device, and basecaller. This information is specified to the tool through a 'model', which is a string of text in the following format:
```
{pore}_{device}_{caller variant}_{caller version}
```

The pipeline requires a medaka model as input. To see a list of medaka models, use the command ```medaka tools list_models```.

As well, it is important to note, the models will *not* contain individual versions of Guppy. Thus, you should choose the version most closest to and less than the version you used (Ex: using Guppy 6.00, enter R***_min_hac_g507)

## Installation
To install this pipeline, enter the following commands:
```
# Clone the repository
git clone https://github.com/rchapman2000/xgen-ont-pipeline

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate xGenONTAssembly
```

## Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd xgen-ont-pipeline

# Use git to pull the latest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate xGenONTAssembly
conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must be either in the same directory as the main.nf file or reference the file location
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --reference REFERENCE_FASTA --primers XGEN_MASTERFILE --model MEDAKA_MODEL
```

### Optional Arguments 
The pipeline also supports the following optional arguments:
| Option | Type | Description |
|---|---|---|
| --minCov | *int* | The minimum depth of coverage required to report a site in the final consensus genome. Sites below this threshold will be masked [Default = 20] |
| --threads | *int* | The number of CPU threads that can be use to run pipeline tools in parallel |

To view the list of options from the command line, use the following command:
```
nextflow run main.nf --help
```