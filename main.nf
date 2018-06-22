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
    .into { samplesGatk ; samplesManta ; samplesStrelka}

process gatk {

	tag { parameters.name }
		     
    input:
    val(parameters) from samplesGatk
    
    output:
    file('*.vcf') into outGatk
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    gatk --spark-runner LOCAL \
    	 --java-options '-Xmx!{task.memory.toGiga()}G' \
    	Mutect2 \
    	-R !{params.ref} \
    	-I !{parameters.tumor} \
    	-I !{parameters.normal} \
    	-tumor !{parameters.name}T \
    	-normal !{parameters.name}N \
    	-O !{parameters.name}.vcf
	
    '''
}

process manta {

	tag { parameters.name }
		     
    input:
    val(parameters) from samplesManta
    
    output:
    file ("manta/results/variants/candidateSmallIndels.vcf.gz*") into outManta
    
    shell:
    '''
        
    shopt -s expand_aliases
    
    configManta.py --normalBam !{parameters.normal} \
    			   --tumorBam !{parameters.tumor} \
    			   --referenceFasta !{params.ref} \
    			   --runDir manta \
    			   --callRegions !{workflow.scriptFile.getParent() + "/util/hg19_chromosomes.bed.gz"}
    			   
    ${PWD}/manta/runWorkflow.py -m local -j !{task.cpus} -g !{task.memory.toGiga()}
	
    '''
}

process strelka {

	tag { parameters.name }
		     
    input:
    val(parameters) from samplesStrelka
    file(mantaVcf) from outManta
    
    output:
    file('strelka/results/variants/*') into outStrelka
    
    shell:
    '''
        
    shopt -s expand_aliases
        
    configureStrelkaSomaticWorkflow.py --normalBam !{parameters.normal} \
    			   --tumorBam !{parameters.tumor} \
    			   --referenceFasta !{params.ref} \
    			   --runDir strelka \
    			   --callRegions !{workflow.scriptFile.getParent() + "/util/hg19_chromosomes.bed.gz"} \
    			   --indelCandidates !{mantaVcf[0]}
    			   
    ${PWD}/strelka/runWorkflow.py -m local -j !{task.cpus} -g !{task.memory.toGiga()}
	
    '''
}
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
