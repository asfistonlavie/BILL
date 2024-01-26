rule sniffles:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_sniffles.vcf"
	log:
		sniffles=f"{{path_res}}/logs/sniffles/{{sample}}.trimed{{trim}}.log",
	benchmark:
		f"{{path_res}}/benchmarks/sniffles/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['sniffles_options']
	shell:
		"""
		sniffles --threads {threads} {params} --reference {input.ref} --sample-id {wildcards.sample}.trimed{wildcards.trim} --input {input.sorted} --vcf {output.vcf} > {log.sniffles} 2>> {log.sniffles}
		"""

rule vcf_sniffles:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_sniffles.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule cutesv:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_cutesv.vcf"
	log:
		cutesv=f"{{path_res}}/logs/cutesv/{{sample}}.trimed{{trim}}.log",
	benchmark:
		f"{{path_res}}/benchmarks/cutesv/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['cutesv_options']
	shell:
		"""
		mkdir -p {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}
		cuteSV --threads {threads} {params} --sample {wildcards.sample}.trimed{wildcards.trim} {input.sorted} {input.ref} {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim} 2>> {log.cutesv}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}
		"""

rule vcf_cutesv:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_cutesv.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule svim:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_svim.vcf"
	log:
		svim=f"{{path_res}}/logs/svim/{{sample}}.trimed{{trim}}.log",
	benchmark:
		f"{{path_res}}/benchmarks/svim/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['svim_options']
	shell:
		"""
		svim alignment {params} --sample {wildcards.sample}.trimed{wildcards.trim} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim {input.sorted} {input.ref} 2> {log.svim}
		bcftools filter --threads 2 --exclude "SUPPORT<10" --output {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim/variants.vcf 2>> {log.svim}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim
		"""

rule vcf_svim:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_svim.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule nanovar:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_nanovar.vcf"
	log:
		nanovar=f"{{path_res}}/logs/nanovar/{{sample}}.trimed{{trim}}.log",
	benchmark:
		f"{{path_res}}/benchmarks/nanovar/{{sample}}.trimed{{trim}}.txt"
	threads: 4
	params:
		config['nanovar_options']
	shell:
		"""
		nanovar --threads {threads} {params} {input.sorted} {input.ref} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar 2> {log.nanovar}
		echo {wildcards.sample}.trimed{wildcards.trim} > {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt
		bcftools reheader --samples {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt --output {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar/{wildcards.sample}.trimed{wildcards.trim}.aligned.sorted.nanovar.pass.vcf 2>> {log.nanovar}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt
		"""

rule vcf_nanovar:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_nanovar.vcf",
			sample=config['samples'],
			trim=config['trim'])

rule nanosv:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		bed=f"{{path_res}}/{{sample}}.trimed{{trim}}.bedgraph"
	output:
		vcf=f"{{path_res}}/{{sample}}.trimed{{trim}}.sv_nanosv.vcf"
	log:
		nanosv=f"{{path_res}}/logs/nanosv/{{sample}}.trimed{{trim}}.log",
	benchmark:
		f"{{path_res}}/benchmarks/nanosv/{{sample}}.trimed{{trim}}.txt"
	threads: 8
	params:
		config['nanosv_options']
	shell:
		"""
		NanoSV --threads {threads} {params} --bed {input.bed} --output {output.vcf} {input.sorted} 2> {log.nanosv}
		bcftools sort --output {output.vcf} {output.vcf} 2>> {log.nanosv} 
		bcftools filter --threads 2 --exclude "SVLEN<10" --output {output.vcf}.tmp {output.vcf} 2>> {log.nanosv}
		bcftools filter --threads 2 --exclude "DV<10" --output {output.vcf} {output.vcf}.tmp 2>> {log.nanosv}
		echo {wildcards.sample}.trimed{wildcards.trim} > {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt
		bcftools reheader --threads 2 --samples {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt --output {output.vcf} {output.vcf}
		rm {output.vcf}.tmp {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt
		"""

rule vcf_nanosv:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.sv_nanosv.vcf",
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
	benchmark:
		f"{{path_res}}/benchmarks/medaka/consensus/{{sample}}.trimed{{trim}}.txt"
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
	benchmark:
		f"{{path_res}}/benchmarks/medaka/snp/{{sample}}.trimed{{trim}}.txt"
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
