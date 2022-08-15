#!/usr/bin/env bash
#BSUB -J TCGA_APA
#BSUB -o log/snakemake_%J.out
#BSUB -e log/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -R "select[hname!=compute03 && hname!=compute13]"

set -o nounset -o pipefail -o errexit -x


args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "select[hname!=compute09 && hname!=compute13 && hname!=compute14] {params.memory}" -n {threads}'

for cancer in read ucec thym cesc pcpg esca paad  chol blca kich stad kirp coad hnsc prad lihc lusc luad thca kirc brca; do
        cat ref/config_template.yaml > config_${cancer}.yaml
	    echo "  \"ref/${cancer}_dict.tsv\"" >> config_${cancer}.yaml
	    echo "" >> config_${cancer}.yaml
	    echo "DIR:" >> config_${cancer}.yaml
	    echo "  \"${cancer^^}\"" >> config_${cancer}.yaml

	cut -f1 ref/${cancer}_dict.tsv | xargs -n10 -d '\n' | tr ' ' ',' | while read i; do
	    echo $(eval echo ${cancer^^}/{$i}\-dapars/{$i}\-dapars_All_Prediction_Results.txt | sed 's/[{}]//g')
	    snakemake --drmaa "$args" --snakefile Snakefile --jobs 40 \
		    --latency-wait 300 --rerun-incomplete --keep-going --configfile config_${cancer}.yaml \
		    --restart-times 2 $(eval echo ${cancer^^}/{$i}-done_dapars.txt | sed 's/[{}]//g')
	done
done
