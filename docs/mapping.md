# Read mapping

<!-- resulmé -->

## Mapping

![mapping_rules](img/mapping.svg "Mapping rules")

**minimap2**

- Input:
	- Fastq || fastq.gz file (reads)  
&emsp;&rarr; provided by `seqkit seq` rule
	- Fasta file (reference)
- Output:
	- SAM file  
&emsp;&rarr; used by `samtools view` rule
- Description:  
&emsp;Mapped the reads on the reference provided.
- Default options:
	- `--MD` &rarr; output the MD tag
	- `-a` &rarr; choose output SAM format
	- `-x map-ont` &rarr; choose Nanopore vs reference mapping

### samtools view

### samtools sort

## Indexing

<!-- samtools index -->

## Statistic control

<!-- seqkit stats -->
<!-- samtools flagstat -->
