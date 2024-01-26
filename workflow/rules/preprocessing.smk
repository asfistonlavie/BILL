rule seqkit_seq:
	input:
		fastq=check_fastq_gz
	output:
		trimed=temp(f"{{path_res}}/{{sample}}.trimed{{trim}}.fastq.gz")
	log:
		f"{{path_res}}/logs/seqkit/seq/{{sample}}.trimed{{trim}}.log"
	benchmark:
		f"{{path_res}}/benchmarks/seqkit/seq/{{sample}}.trimed{{trim}}.txt"
	threads: 2
	params:
		config['seqkit_seq_options']
	shell:
		"seqkit seq -j {threads} {params} -m {wildcards.trim} {input} -o {output} 2> {log}"
    
rule trim:
	input:
		expand(
			f"{config['res_dir']}/{{sample}}.trimed{{trim}}.fastq.gz",
			sample=config['samples'],
			trim=config['trim']
		)
