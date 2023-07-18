rule sniffles:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv.vcf"
	log:
		sniffles=f"{{path_res}}/logs/sniffles/{{sample}}.trimed{{trim}}.log",
		bcftools=f"{{path_res}}/logs/bcftools/reheader/{{sample}}.trimed{{trim}}.log"
	threads: 4
	params:
		config['sniffles_options']
	shell:
		"""
		sniffles -t {threads} {params} -i {input.sorted} -v {output}.tmp > {log.sniffles} 2>> {log.sniffles}
		echo {wildcards.sample}.trimed{wildcards.trim} > {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name.txt
		bcftools reheader -s {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name.txt -o {output} {output}.tmp 2> {log.bcftools}
		rm {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name.txt {output}.tmp
		"""

rule vcf:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule medaka_consensus:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output
	output:
		hdf=f"{{path_res}}/{{sample}}.trimed{{trim}}.hdf"
	log:
		f"{{path_res}}/logs/medaka/consensus/{{sample}}.trimed{{trim}}.log"
	threads: 2
	params:
		config['medaka_consensus_options']
	shell:
		"medaka consensus --threads {threads} {params} {input.sorted} {output} 2> {log}"

rule hdf:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.hdf",
			sample=config['samples'],
			trim=config['trim'])

rule medaka_snp:
	input:
		ref=f"{config['ref_dir']}/{config['ref_name']}",
		hdf=rules.medaka_consensus.output
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.snp.vcf"
	log:
		f"{{path_res}}/logs/medaka/snp/{{sample}}.trimed{{trim}}.log"
	threads: 4
	params:
		config['medaka_snp_options']
	shell:
		"medaka snp {params} {input} {output} 2> {log}"

rule snp:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.snp.vcf",
			sample=config['samples'],
			trim=config['trim'])
