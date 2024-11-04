# BioInformatics Learning Lab - BILL

<img src="img/logo_bill.jpg" alt="Logo BILL" title="Logo BILL" align="left" width=250 height=250/>

## Variant calling pipeline designed for teaching

The Bioinformatics Learning Lab (BILL) is a teaching unit of the Master of Bioinformatics of the University of Montpellier. Students take part in a research project analysing structural variants (SVs) and small nucleotide variants (SNVs).

They perform DNA extraction, sequencing, data analysis and interpretation of the results. 

This pipeline starts by trimming the read files by removing reads smaller than 1,000 bp. It then proceeds to align the reads against the virus genomic reference. It processes the alignment by removing unaligned reads and converting them to a sorted binary format. It then performs a variant calling step and filters the resulting variants. Finally, it merges all variant files into a _VCF_ file. Some statistical commands appear throughout the pipeline to check the quality of the data or results.

<img src="img/pipeline.svg" alt="Pipeline" title="Pipeline"/>

- [Getting Started](#getting-started)
    - [Dependencies](#dependencies)
    - [How to install it](#install)
    - [How to use it](#how-to-use-it)
- [Manual](https://souliera.github.io/BILL/)
- [Authors](#authors)
- [Acknowledgements](#acknowledgement)

## Getting Started

### Dependencies

- snakemake v7.21.0
- seqkit v2.3.0
- minimap2 v2.24-r1122
- samtools v1.16.1
- bamCoverage v3.5.1
- plotCoverage v3.5.1
- sniffles2 v2.2
- tabix v1.16
- bgzip v1.16
- bcftools v1.16
- medaka v1.11.3

### How to install it

Clone the repository wherever you want on your local:

```
git clone git@gitlab.com:asfistonlavie/bill.git # by SSH 
git clone https://gitlab.com/asfistonlavie/bill.git # by HTTPS
```

### How to use it

There are two ways to run this pipeline on a cluster (i.e NGSTC): with the srun command (`srun snakemake --cores <nb_core_max> --configfile workflow/config/config.yaml [options] <target>`), or with a cluster config (`snakemake --profile workflow/ngstc --profile workflow/config/config.yaml [options] <target>`).

There are three targets to use this pipeline: (1) a file name, (2) a rule name, or (3) the whole pipeline.

If you just want a specific file, run:
```
srun snakemake --cores <nb_core_max> --configfile workflow/config/config.yaml -k -p <file_name>
```
It will automatically find the correct rule to run based on the file name. File names are constrained by the snakefile (see the [wiki](https://gitlab.com/souliera/bill/-/wikis/Release-2024/Rule-details) for correct file name format).

If you want all of a type of file, run :
```
srun snakemake --cores <nb_core_max> --configfile workflow/config/config.yaml -k -p <rule_name>
```

## How to cite the pipeline

Variant calling pipeline designed for the Bioinformatics Learning Lab (BILL: release_2023.1)
A Soulier, A Arnoux, C Breton, E Cherif, AS Fiston-Lavier
Zenodo. https://doi.org/10.5281/zenodo.10020027

## Authors

- Arnaud Soulier [Github profil](https://github.com/souliera)
- Catherine Breton [Github profil](https://github.com/CathyBreton)
- Emira Cherif [Github profil](https://github.com/emiracherif)
- Anna-Sophie Fiston-Lavier [Github profil](https://github.com/asfistonlavie)


## Acknowledgement
 
We would like to thank the people involved in the BILL project, Jean-Christophe Avarre, Anne-Sophie Gosselin-Grenet and Marie-Ka Tilak, as well as all the students who took part in the project. 

[![Website](https://img.shields.io/website?up_message=up&up_color=green&down_message=down&down_color=red&url=https://informatique-fds.edu.umontpellier.fr/etudiants/masters-transdisciplinaires/master-bioinformatique/bill-bioinformatics-learning-lab/)](https://informatique-fds.edu.umontpellier.fr/etudiants/masters-transdisciplinaires/master-bioinformatique/bill-bioinformatics-learning-lab/)
![GitHub](https://img.shields.io/github/license/asfistonlavie/BILL)
![GitHub contributors](https://img.shields.io/github/contributors/asfistonlavie/BILL)
[![Github release](https://badgen.net/github/releases/asfistonlavie/BILL)](https://github.com/asfistonlavie/BILL/releases)
[![Static Badge](https://img.shields.io/badge/wiki-yes-green)](https://github.com/asfistonlavie/BILL/wiki)
