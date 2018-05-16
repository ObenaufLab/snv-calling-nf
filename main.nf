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
        					
        					The following columsn are required:
        						- name: Name of the sample
        						- normal: Path to normal bam file
        						- tumor: Path to tumor bam file
        						
        --ref				Path to reference fasta file

    Profiles:
        standard            local execution
        sge			        SGE execution with singularity on IMPIMBA1
        ii2                 SLURM execution with singularity on IMPIMBA2
        
    Docker:
    obenauflab/variant-circos-nf:latest
    
    Author:
    Tobias Neumann (tobias.neumann@imp.ac.at)
    """.stripIndent()
}

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .set { samples }

process gatk {

	tag { name }
     
    input:
    val(parameters) from samples
    
    shell:
    ''' 
    shopt -s expand_aliases
	
	echo !{parameters.name}
	echo !{parameters.normal}
	echo !{parameters.tumor}
	
    '''
}
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
