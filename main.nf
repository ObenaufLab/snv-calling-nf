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

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .into { samplesChannel ; setupChannel }

process somaticSeqSetup {

	tag { parameters.name }
		     
    input:
    val(parameters) from samplesChannel
    
    output:
    file('somaticseq') into outSomaticSeqSetup
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    /opt/somaticseq/utilities/dockered_pipelines/submit_callers_multiThreads.sh \
    	-o somaticseq \
    	--normal-bam !{parameters.normal} \
    	--tumor-bam !{parameters.tumor} \
    	--human-reference !{params.ref} \
    	--action echo --mutect2 --somaticsniper --vardict --scalpel --strelka --somaticseq --lofreq \
    	--threads !{task.cpus} --dbsnp /groups/zuber/zubarchive/USERS/tobias/hg38/GATK/Homo_sapiens_assembly38.dbsnp138.vcf
	
    '''
}

process lofreq {

	tag { parameters.name }
		     
    input:
    val(parameters) from samplesChannel
    file(somaticseq) from outSomaticSeqSetup
    
    output:
    file('somaticseq') into out
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    sed -i '/^#/d' somaticseq/1/logs/lofreq_*.cmd
    sed -i 's/\\/mnt\\///g' somaticseq/1/logs/lofreq_*.cmd
    sed -i 's/.*lethalfang/singularity exec \\/groups\\/zuber\\/zubarchive\\/USERS\\/tobias\\/.singularity/g' somaticseq/1/logs/lofreq_*.cmd
    sed -i 's/\\.singularity\\S*/&.img/' somaticseq/1/logs/lofreq_*.cmd
    echo -e "#SBATCH --mem 49152\n$(cat somaticseq/1/logs/lofreq_*.cmd)"
    echo -e "#SBATCH --error=somaticseq/1/logs/lofreq_slurm-%j.err\n$(cat somaticseq/1/logs/lofreq_*.cmd)"
    echo -e "#SBATCH --output=somaticseq/1/logs/lofreq_slurm-%j.out\n$(cat somaticseq/1/logs/lofreq_*.cmd)"
    echo -e "#!/usr/bin/env bash\\n$(cat somaticseq/1/logs/lofreq_*.cmd)"
	
    '''
}

 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
