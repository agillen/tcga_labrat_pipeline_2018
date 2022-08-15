shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; ")
import glob
import os

''' rules to process TCGA APA data '''

# configfile: "config.yaml"
  
KEY = config["KEY"]
SIZES = config["SIZES"]
TCGA_DICT = config["TCGA_DICT"]
SALMON = config["SALMON"]
SALMON_IDX = config["SALMON_IDX"]
DAPARS_IDX = config["DAPARS_IDX"]
DAPARS = config["DAPARS"]
DAPARS_CFG = config["DAPARS_CFG"]
DIR = config["DIR"]

d = {}
with open(TCGA_DICT) as f:
    for line in f:
       (key, id_norm, id_tumor, file_norm, file_tumor) = line.split()
       d[key] = [id_norm, id_tumor, file_norm, file_tumor]

rule all:
	input: expand('{DIR}/{patient}-done.txt', patient = [*d], DIR=DIR)

rule gdc_client:
	output: "{DIR}/{patient}-{sample}_pass.out"
	params:
	    job_name = '{patient}-{sample}_.gdc_client',
	    memory = "select[mem>5] rusage[mem=5] span[hosts=1]"
	log:
	    'log/{patient}-{sample}_gdc-client'
	run:
	    if wildcards.sample == 'norm':
	    	input_gdc = d[wildcards.patient][0]
	    else:
	    	input_gdc = d[wildcards.patient][1]
	    shell("module load python/2.7.8")
	    shell("rm -rf {DIR}/{input_gdc}")
	    shell("~/build/gdc-client/bin/gdc-client download {input_gdc} -t {KEY} -d {DIR}")
	    shell("touch {output}")
	    
rule bedgraph:
	input: rules.gdc_client.output
	output: "{DIR}/{patient}-{sample}.bedgraph"
	params:
	    job_name = '{patient}-{sample}_bedgraph',
	    memory = "select[mem>5] rusage[mem=10] span[hosts=1]"
	log:
	    'log/{patient}-{sample}_bedgraph'
	run:
	    if wildcards.sample == 'norm':
	    	input_bam = DIR+'/'+d[wildcards.patient][0]+'/'+d[wildcards.patient][2]
	    else:
	    	input_bam = DIR+'/'+d[wildcards.patient][1]+'/'+d[wildcards.patient][3]
	    shell("genomeCoverageBed -bg -ibam {input_bam} -g {SIZES} -split > {output}")

rule sort_bam:
        input: rules.gdc_client.output
        output: "{DIR}/{patient}-{sample}_sorted.bam"
        params:
            sample = "{DIR}/{patient}-{sample}_sorted",
            job_name = '{patient}-{sample}_sort_Bed',
            memory = "select[mem>34] rusage[mem=34] span[hosts=1]"
        log:
            'log/{patient}-{sample}_sort_bed'
        threads: 8
        run:
            if wildcards.sample == 'norm':
                input_bam = DIR+'/'+d[wildcards.patient][0]+'/'+d[wildcards.patient][2]
            else:
                input_bam = DIR+'/'+d[wildcards.patient][1]+'/'+d[wildcards.patient][3]
            shell("samtools sort -n -@ {threads} -m 4G {input_bam} -T {params.sample} > {output}")

rule bamtofastq:
	input: "{DIR}/{patient}-{sample}_sorted.bam"
	output:
	    read0 = "{DIR}/{patient}-{sample}_sorted_r0.fastq.gz",
	    read1 = "{DIR}/{patient}-{sample}_sorted_r1.fastq.gz",
	    read2 = "{DIR}/{patient}-{sample}_sorted_r2.fastq.gz",
	    paired = "{DIR}/{patient}-{sample}_paired.txt"
	params:
	    job_name = '{patient}-{sample}_bamtofastq',
	    memory = "select[mem>20] rusage[mem=20] span[hosts=1]"
	log:
	    'log/{patient}-{sample}_bamtofastq'
	shell:
	    """
	    samtools fastq -0 {output.read0} -1 {output.read1} -2 {output.read2} {input}
	    if [ -z $(zcat {output.read1} | head -c1) ]; then
	        touch {output.paired}
	    else
	        echo 1 > {output.paired}
	    fi
        rm -f {input}
	"""

rule salmon:
	input:
	    read0 = "{DIR}/{patient}-{sample}_sorted_r0.fastq.gz",
	    read1 = "{DIR}/{patient}-{sample}_sorted_r1.fastq.gz",
	    read2 = "{DIR}/{patient}-{sample}_sorted_r2.fastq.gz",
	    paired_norm = "{DIR}/{patient}-norm_paired.txt",
	    paired_tumor = "{DIR}/{patient}-tumor_paired.txt"
	output: "{DIR}/{patient}-{sample}_salmon/quant.sf"
	params:
	    job_name = '{patient}-{sample}_salmon',
	    memory = "select[mem>5] rusage[mem=5] span[hosts=1]",
	    sample = "{DIR}/{patient}-{sample}_salmon",
	    unpaired = "{DIR}/{patient}-unpaired.txt",
	    paired = "{DIR}/{patient}-{sample}_paired.txt"
	log:
	    'log/{patient}-{sample}_salmon'
	threads: 8
	shell:
	    """
	    if [ -s {input.paired_norm} ] && [ -s {input.paired_tumor} ]; then
	        {SALMON} quant --libType A --seqBias --gcBias \
	        -1 {input.read1} -2 {input.read2} \
	        -o {params.sample} -p 8 --index {SALMON_IDX}
	    else
	        if [ -s {params.paired} ] ; then
	            read0={input.read1}
	        else
	            read0={input.read0}
	        fi
	        {SALMON} quant --libType A --seqBias --gcBias \
	        -r ${{read0}} \
	        -o {params.sample} -p 8 --index {SALMON_IDX}
	        touch {params.unpaired}
	    fi
	    rm -f {input.read0} {input.read1} {input.read2}
        """

rule dapars:
	input:
	    norm = "{DIR}/{patient}-norm.bedgraph",
	    tumor = "{DIR}/{patient}-tumor.bedgraph"
	output: "{DIR}/{patient}-dapars/{patient}-dapars_All_Prediction_Results.txt"
	params:
	    job_name = '{patient}-dapars',
	    memory = "select[mem>5] rusage[mem=5] span[hosts=1] -c 240",
	    config = "{DIR}/{patient}-dapars.config",
	    out = "{patient}-dapars"
	log:
	    'log/{patient}-dapars'
	shell:
	    """
	    module load python
	    echo Annotated_3UTR={DAPARS_IDX} > {params.config}
	    echo Group1_Tophat_aligned_Wig={input.norm} >> {params.config}
	    echo Group2_Tophat_aligned_Wig={input.tumor} >> {params.config}
	    echo Output_directory={DIR}/{params.out}/ >> {params.config}
	    echo Output_result_file={params.out} >> {params.config}
	    cat {DAPARS_CFG} >> {params.config}
	    python {DAPARS} {params.config}
	    rm -f {input.norm} {input.tumor} {params.config}
	"""

rule dummy:
	input:
#	    expand("{{DIR}}/{{patient}}-{sample}_salmon/quant.sf", sample = ['norm', 'tumor']),
	    rules.dapars.output
#            norm = "{DIR}/{patient}-norm.bedgraph",
#            tumor = "{DIR}/{patient}-tumor.bedgraph"
	output: "{DIR}/{patient}-done_dapars.txt"
	params:
	    job_name = '{patient}-dummy',
	    memory = "select[mem>5] rusage[mem=5] span[hosts=1]",
	    paired_txt = "{DIR}/{patient}-*_paired.txt",
	    pass_txt = "{DIR}/{patient}-*pass.out"
	log:
	    'log/{patient}-dummy'
	run:
	    input_norm = DIR+'/'+d[wildcards.patient][0]
	    input_tumor = DIR+'/'+d[wildcards.patient][1]
	    shell("rm -rf {params.pass_txt} {input_norm} {input_tumor} {params.paired_txt}") 
	    shell("touch {output}")
