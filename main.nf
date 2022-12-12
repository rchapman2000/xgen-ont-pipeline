#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2
def helpMessage() {
    log.info"""
Referenced-Based ONT xGen Amplicon Sequencing Assembly Pipeline

This pipeline was built specifically for the adaptation of IDT's xGen
Amplicon panels to Oxford Nanopore Sequencing.
It takes a set of fastq files produced by ONT Sequencing, and generates a
consensus genome from the data. The pipeline begins by aligning reads using
minimap2 and then soft-clips primer sequences using the tool Primerclip which is designed
for xGen technology. Next, the alignment is corrected using medaka, and variants are called
from this correction using Medaka and Longshot. Finally, depth masking and variants are applied
to the reference genome.

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --reference REFERENCE_FASTA --primers XGEN_MASTERFILE --model MEDAKA_MODEL

OPTIONS:

--input INPUT_DIR - [Required] A directory containing nanopore fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--reference REFERENCE_FASTA - [Required] A reference genome to align reads to.

--primers XGEN_MASTERFILE - [Required] The masterfile supplied by IDT specifying primer locations to be clipped using PrimerClip

--model MEDAKA_MODEL - [Required] Medaka requires a 'model' which corresponds to the flow-cell type/basecaller parameters used to correct for errors typical to that technology. Use 'medaka tools list_models' to find this.

OPTIONAL:
    
    --minCov INT - The minimum coverage below which a position will be masked [Default = 20]

    --threads INT - The number of CPU threads that can be use to run pipeline tools in parallel
"""
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.reference = false
params.output = false
params.primers = false
params.threads = 1
params.minCov = 20
params.model = false

println "Input Directory: ${params.input}"

include { Setup } from './modules.nf'
include { MiniMap2_Alignment } from './modules.nf'
include { XGen_Primer_Clip } from './modules.nf'
include { Medaka_Consensus } from './modules.nf'
include { Call_Variants } from './modules.nf'
include { Generate_Consensus_Genome } from './modules.nf'
include { Write_Summary } from './modules.nf'

// Checks the input parameter
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromPath("${params.input}*.fastq*")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it.getSimpleName(), it]}


// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
    println(outDir)
}

// Checks the reference parameter. For this, we cannot use an
// input channel like was used for the input files. Using an input channel
// will cause Nextflow to only iterate once as the reference 
// channel would only only have 1 file in it. Thus, we manually parse
// the reference file into a tuple.
refData = ''
refName = ''
if (params.reference == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: no reference file proivded. Pipeline requires a reference file."
    exit(1)
}
else if (!(file(params.reference).exists())) {
    // If the reference file provided does not exist, notify the user and exit.
    println "ERROR: ${params.reference} does not exist."
    exit(1)
}
else {
    // Process the reference file to be supplied to the index step.
    
    // Parse the file provided into a file object.
    ref = file(params.reference)

    // Grab the basename of the file.
    refName = ref.getBaseName()

    // Place the file basename and file object into
    // a tuple.
    refData = tuple(refName, ref)
}
println(refData)

// Process the xGen primers parameter.
primerfile = ''
primerFileName = ''
if (params.primers != false) {
    // If the parameter is provided, check if the file exists.
    if (file(params.primers).isFile())
    {
        // If the file does exist, parse it into
        // a file object.
        primerfile = file(params.primers)
        primerFileName = primerfile.getName()
    }
    else {
        // If the file does not exist, notify the user and exit.
        println "ERROR: ${params.primers} does not exist."
        exit(1)
    }
}
else {
    // If no primer file was supplied, notify the user and exit
    println "ERROR: no primer masterfile provided. Pipeline requires a primer file."
    exit(1)
}

// Processes the medaka model parameter.
model = ''
if (params.model == false) {
    // If no model was provided, notify the user, provide them instructions,
    // on how to find this value from medaka, and exit.
    println "ERROR: no ONT model provided. Medaka requires a model in the format:"
    println ""
    println "{pore}_{device}_{caller variant}_{caller version}"
    println ""
    println "To see existing models enter: medaka tools list_models"
    exit(1)
}
else {
    // If the parameter was provided. Store it in a variable.
    model = params.model
}

// To keep track of the summary, a string will be passed between each process that contains
// metrics generated by each. At the end of each process, the metrics collected will be
// appended to the end and returned as output from that process. At the end of the workflow,
// the metrics will be written to the summary file.
workflow {
    // Creates files to store analysis parameters and statistics.
    Setup ( refName, params.minCov, primerFileName, model, outDir )

    // Maps nanopore reads to the provided reference using medaka.
    MiniMap2_Alignment( inputFiles_ch, outDir, refData )

    // Clips primers from the alignment using PrimerCip
    XGen_Primer_Clip( MiniMap2_Alignment.out[0], primerfile, outDir, params.threads, MiniMap2_Alignment.out[1] )

    // Polishes the alignment using medaka consensus.
    Medaka_Consensus( XGen_Primer_Clip.out[0], model, outDir, XGen_Primer_Clip.out[1] )

    // Calls variants using medaka and longshot and filters them.
    Call_Variants( Medaka_Consensus.out[0], baseDir, outDir, refData, params.minCov, Medaka_Consensus.out[1] )

    // Generates a consensus genome by performing a depth mask and then applying variants to
    // the reference.
    Generate_Consensus_Genome( Call_Variants.out[0], baseDir, outDir, refData, params.minCov, Call_Variants.out[2] )

    // Writes statistics to an analysis summary file.
    Write_Summary( Generate_Consensus_Genome.out[2], outDir )
}