##############
### INPUTS ###
##############

rule seqkit_stats_raw:
	input:
		fastq=check_fastq_gz
	output:
		stats=protected(f"{{path_res}}/{{sample}}.seqstats")
	log:
		f"{{path_res}}/logs/seqkit/stats/{{sample}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/seqkit/stats/{{sample}}.txt"
	threads: 2
	params:
		config['seqkit_stats_options']
	shell:
		"seqkit stats -j {threads} {params} {input} -o {output} 2> {log}"

rule seq_stats_raw:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.seqstats",
			sample=config['samples']
		)

rule seqkit_stats_trimed:
	input:
		fastq=rules.seqkit_seq.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.seqstats")
	log:
		f"{{path_res}}/logs/seqkit/stats/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/seqkit/stats/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['seqkit_stats_options']
	shell:
		"seqkit stats -j {threads} {params} {input} -o {output} 2> {log}"

rule seq_stats_trimed:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.seqstats",
			sample=config['samples'],
			trim=config['trim']
		)


#################
### ALIGNMENT ###
#################

rule samtools_flagstat:
	input:
		sam=rules.minimap2.output
	output:
		flagstat=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.flagstat")
	log:
		f"{{path_res}}/logs/samtools/flagstat/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/samtools/flagstat/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['samtools_flagstat_options']
	shell:
		"samtools flagstat -@ {threads} {params} {input} > {output} 2> {log}"

rule flagstat:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.flagstat",
			sample=config['samples'],
			trim=config['trim'])

rule plot_flagstat:
	input:
		flagstat=get_flagstat,
		script=f"workflow/scripts/plot_flagstat.py"
	output:
		png=f"{{path_res}}/{{group}}.trimed{{trim}}.flagstat.png"
	log:
		f"{{path_res}}/logs/python/flagstat/{{group}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/python/flagstat/{{group}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['python_flagstat_plot_options']
	shell:
		"python {input.script} {input.flagstat} {output.png}"

rule flagstat_plot:
	input:
		expand(
			f"{config['res_dir']}/{{group}}.trimed{{trim}}.flagstat.png",
			group=config['group'],
			trim=config['trim']
		)

########################
### DEPTH & COVERAGE ###
########################

rule bamCoverage:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		bedgraph=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.bedgraph")
	log:
		f"{{path_res}}/logs/deeptools/bamCoverage/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/bamCoverage/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['bamCoverage_options']
	shell:
		"bamCoverage -p {threads} --effectiveGenomeSize `grep -v \"^>\" {input.ref} | tr -d \"\r\n\" | wc -c` {params} -b {input.sorted} -o {output} 2> {log}"

rule bedgraph:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.bedgraph",
			sample=config['samples'],
			trim=config['trim'])

rule plotCoverage:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output
	output:
		plot=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.pdf")
	log:
		f"{{path_res}}/logs/deeptools/plotCoverage/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/plotCoverage/{{sample}}.trimed{{trim}}.txt"
	threads: 8
	params:
		config['plotCoverage_options']
	shell:
		"plotCoverage -p {threads} {params} -b {input.sorted} -o {output} > {log} 2>> {log}"


rule plot:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.pdf",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_sv_sniffles:
	input:
		vcf=rules.sniffles.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_sniffles.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.sv_sniffles.log"
	benchmark:
		f"{{path_res}}/benchmarks/bcftools/stats/{{sample}}.trimed{{trim}}.sv_sniffles.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule sv_stats_sniffles:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_sniffles.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_sv_cutesv:
	input:
		vcf=rules.cutesv.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_cutesv.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.sv_cutesv.log"
	benchmark:
		f"{{path_res}}/benchmarks/bcftools/stats/{{sample}}.trimed{{trim}}.sv_cutesv.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule sv_stats_cutesv:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_cutesv.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_sv_nanosv:
	input:
		vcf=rules.nanosv.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_nanosv.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.sv_nanosv.log"
	benchmark:
		f"{{path_res}}/benchmarks/bcftools/stats/{{sample}}.trimed{{trim}}.sv_nanosv.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule sv_stats_nanosv:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_nanosv.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_sv_svim:
	input:
		vcf=rules.svim.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_svim.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.sv_svim.log"
	benchmark:
		f"{{path_res}}/benchmarks/bcftools/stats/{{sample}}.trimed{{trim}}.sv_svim.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule sv_stats_svim:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_svim.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_sv_nanovar:
	input:
		vcf=rules.nanovar.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_nanovar.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.sv_nanovar.log"
	benchmark:
		f"{{path_res}}/benchmarks/nanovar/stats/{{sample}}.trimed{{trim}}.sv_nanovar.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule sv_stats_nanovar:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_nanovar.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_stats_snp:
	input:
		vcf=rules.medaka_snp.output
	output:
		stats=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.snp.vcf.stats")
	log:
		f"{{path_res}}/logs/bcftools/stats/{{sample}}.trimed{{trim}}.snp.log"
	benchmark:
		f"{{path_res}}/benchmarks/bcftools/stats/{{sample}}.trimed{{trim}}.snp.txt"
	threads: 2
	params:
		config['bcftools_stats_options']
	shell:
		"bcftools stats --threads {threads} {params} {input} > {output} 2> {log}"

rule snp_stats:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.snp.vcf.stats",
			sample=config['samples'],
			trim=config['trim'])
