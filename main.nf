#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2018 Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

def helpMessage() {
    log.info"""
    ================================================================
    snv-calling-nf
    ================================================================
    DESCRIPTION
    Usage:
    nextflow run obenauflab/snv-calling-nf
    Options:
        --samples        	Tab-delimited text file specifying the samples to be 
        					analyzed.
        					
        					The following columns are required:
        						- name: Name of the sample
        						- normal: Path to normal bam file
        						- tumor: Path to tumor bam file
        						
        --ref				Path to reference fasta file

    Profiles:
        standard            local execution
        sge			        SGE execution with singularity on IMPIMBA1
        ii2                 SLURM execution with singularity on IMPIMBA2
        
    Docker:
    broadinstitute/gatk:4.0.4.0
    obenauflab/strelka:latest
    
    Author:
    Tobias Neumann (tobias.neumann@imp.ac.at)
    """.stripIndent()
}

params.submit = false 

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .set { setupChannel }

process somaticSeqSetup {

	tag { parameters.name }
		     
    input:
    val(parameters) from setupChannel
    
    output:
    file("somaticseq_${parameters.name}") into outSomaticSeqSetup
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    /opt/somaticseq/utilities/dockered_pipelines/submit_callers_multiThreads.sh \
    	-o somaticseq_!{parameters.name} \
    	--normal-bam !{parameters.normal} \
    	--tumor-bam !{parameters.tumor} \
    	--human-reference !{params.ref} \
    	--action echo --mutect2 --somaticsniper --vardict --scalpel --strelka --somaticseq --lofreq \
    	--threads !{task.cpus} --dbsnp /groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sapiens_assembly38.dbsnp138.vcf
	
    for log in `ls somaticseq_!{parameters.name}/*/logs/*lofreq*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --qos=long\n$(cat $log)" > $log
    	echo -e "#SBATCH --mem 49152\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/lofreq_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/lofreq_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/*/logs/*mutect2*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/`.*lethalfang/`singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/.*broadinstitute/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH -c 24\n$(cat $log)" > $log
    	echo -e "#SBATCH --mem 49152\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/mutect2_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/mutect2_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/*/logs/*scalpel*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --qos=long\n$(cat $log)" > $log
    	echo -e "#SBATCH --mem 49152\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/scalpel_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/scalpel_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/*/logs/*strelka*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --mem 6144\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/strelka_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/strelka_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/*/logs/*vardict*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --mem 6144\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/vardict_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/vardict_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/logs/*somaticsniper*.cmd`
	do
		filepath=$(dirname $log)
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --mem 6144\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/somaticsniper_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/somaticsniper_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	for log in `ls somaticseq_!{parameters.name}/*/SomaticSeq/logs/sseq_*.cmd`
	do
		filepath=$(dirname $log)
		sed -i '/docker pull/d' $log
    	sed -i '/^#/d' $log
    	sed -i 's/\\/mnt\\///g' $log
    	sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' $log
    	sed -i 's/\\.singularity\\S*/&.img/' $log
    	echo -e "#SBATCH --mem 6144\n$(cat $log)" > $log
    	echo -e "#SBATCH --error=$filepath/sseq_slurm-%j.err\n$(cat $log)" > $log
    	echo -e "#SBATCH --output=$filepath/sseq_slurm-%j.out\n$(cat $log)" > $log
    	echo -e "#!/usr/bin/env bash\\n$(cat $log)" > $log
	done
	
	echo "cd .." >> somaticseq_!{parameters.name}/submitSomaticSeq.sh
	echo "for log in \\`ls somaticseq_!{parameters.name}/*/SomaticSeq/logs/sseq_*.cmd\\`" >> somaticseq_!{parameters.name}/submitSomaticSeq.sh
	echo "do" >> somaticseq_!{parameters.name}/submitSomaticSeq.sh
	echo "	sbatch \\$log" >> somaticseq_!{parameters.name}/submitSomaticSeq.sh
	echo "done" >> somaticseq_!{parameters.name}/submitSomaticSeq.sh
    

    '''
}

 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
