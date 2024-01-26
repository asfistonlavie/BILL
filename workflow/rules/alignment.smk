################
### MINIMAP2 ###
################

rule minimap2:
	input:
		reads=rules.seqkit_seq.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		sam=temp(f"{{path_res}}/{{sample}}.trimed{{trim}}.sam")
	log:
		f"{{path_res}}/logs/minimap/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/minimap/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['minimap2_options']
	shell:
		"minimap2 -t {threads} {params} {input.ref} {input.reads} -o {output} 2> {log}"

rule sam:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sam",
			sample=config['samples'],
			trim=config['trim'])

################
### SAMTOOLS ###
################

rule samtools_view:
	input:
		sam=rules.minimap2.output
	output:
		aligned=temp(f"{{path_res}}/{{sample}}.trimed{{trim}}.aligned.bam")
	log:
		f"{{path_res}}/logs/samtools/view/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/samtools/view/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['samtools_view_options']
	shell:
		"samtools view -@ {threads} {params} {input} -o {output} 2> {log}"

rule aligned:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.aligned.bam",
			sample=config['samples'],
			trim=config['trim'])

rule samtools_sort:
	input:
		aligned=rules.samtools_view.output
	output:
		sorted=f"{{path_res}}/{{sample}}.trimed{{trim}}.aligned.sorted.bam"
	log:
		f"{{path_res}}/logs/samtools/sort/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/samtools/sort/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['samtools_sort_options']
	shell:
		"samtools sort -@ {threads} {params} {input} -o {output} 2> {log}"

rule sorted:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.aligned.sorted.bam",
			sample=config['samples'],
			trim=config['trim'])

rule samtools_index:
	input:
		sorted=rules.samtools_sort.output
	output:
		index=f"{{path_res}}/{{sample}}.trimed{{trim}}.aligned.sorted.bam.bai"
	log:
		f"{{path_res}}/logs/samtools/index/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/samtools/index/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['samtools_index_options']
	shell:
		"samtools index -@ {threads} {params} {input} 2> {log}"

rule align_index:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.aligned.sorted.bam.bai",
			sample=config['samples'],
			trim=config['trim'])

include: "variants_calling.smk"

rule medaka_stitch:
	input:
		hdf=rules.medaka_consensus.output,
		draft=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		fasta=f"{{path_res}}/{{sample}}.trimed{{trim}}.fasta"
	log:
		f"{{path_res}}/logs/medaka/stitch/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/medaka/stitch/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['medaka_stitch_options'],
	shell:
		"medaka stitch --threads {threads} {params} {input.hdf} {input.draft} {output} 2> {log}"

rule assembly:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.fasta",
			sample=config['samples'],
			trim=config['trim'])
