rule bcftools_filter:
	input:
		vcf=rules.sniffles.output
	output:
		filter=protected(f"{{path_res}}/{{sample}}.trimed{{trim}}.sv.filtered.vcf")
	log:
		f"{{path_res}}/logs/bcftools/filter/{{sample}}.trimed{{trim}}.log"
	threads: 2
	params:
		config['bcftools_filter_options']
	shell:
		"bcftools filter --threads {threads} {params} {input} -o {output} 2> {log}"

rule vcf_filter:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv.filtered.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule bgzip:
	input:
		vcf=rules.bcftools_filter.output
	output:
		gz=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv.filtered.vcf.gz"
	log:
		f"{{path_res}}/logs/bgzip/{{sample}}.trimed{{trim}}.log"
	threads: 2
	params:
		config['bgzip_options']
	shell:
		"bgzip -@ {threads} {input} 2> {log}"

rule vcf_gz:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv.filtered.vcf.gz",
			sample=config['samples'],
			trim=config['trim'])

rule tabix:
	input:
		vcf_gz=rules.bgzip.output
	output:
		index=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv.filtered.vcf.gz.tbi"
	log:
		f"{{path_res}}/logs/tabix/{{sample}}.trimed{{trim}}.log"
	threads: 2
	params:
		config['tabix_options']
	shell:
		"tabix {input} 2> {log}"

rule vcf_index:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv.filtered.vcf.gz.tbi",
			sample=config['samples'],
			trim=config['trim'])

rule bcftools_merge:
	input:
		vcfs=get_vcfs,
		index=get_vcfs_index
	output:
		merge=f"{{path_res}}/{{group}}.trimed{{trim}}.sv.filtered.merge.vcf"
	log:
		f"{{path_res}}/logs/bcftools/merge/{{group}}.trimed{{trim}}.log"
	threads: 2
	params:
		config['bcftools_merge_options']
	shell:
		"bcftools merge --threads {threads} {params} {input.vcfs} -o {output} 2> {log}"

rule vcf_merge:
	input:
		expand(
			f"{config['res_dir']}/{{group}}.trimed{{trim}}.sv.filtered.merge.vcf",
			group=config['group'],
			trim=config['trim'])
