def check_fastq_gz(wildcards):
	if Path(f"{config['input_dir']}/{wildcards.sample}.fastq.gz").is_file():
		return f"{config['input_dir']}/{wildcards.sample}.fastq.gz"
	else:
		return f"{config['input_dir']}/{wildcards.sample}.fastq"

def get_flagstat(wildcards):
	return expand(
		f"{wildcards.path_res}/{{sample}}.trimed{wildcards.trim}.flagstat",
		sample=config['samples']
	)

def get_vcfs(wildcards):
	return expand(
		f"{wildcards.path_res}/{{sample}}.trimed{wildcards.trim}.sv.filtered.vcf.gz",
		sample=config['samples']
	)

def get_vcfs_index(wildcards):
	return expand(
		f"{wildcards.path_res}/{{sample}}.trimed{wildcards.trim}.sv.filtered.vcf.gz.tbi",
		sample=config['samples']
	)
