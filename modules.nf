// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        // The name of the reference supplied (for use
        // in the parameters file).
        val refName
        // The minimum coverage below which a site is masked (for
        // use in the parameters file).
        val minCov
        // The name of the primer file supplied.
        val primerFileName
        // The provided medaka model.
        val model
        // The output directory to be used.
        val outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. The name of the reference supplied.
        2. The minimum coverage threshold used for masking.
        3. The name of the primer file provided.
        4. The medaka model provided.

    The summary file will always contain:
        1. Sample name
        2. Raw Reads
        3. Mapped Reads
        4. Clipped, Mapped Reads
        5. Average Read Depth
        6. SNPs Passing Filtering
        7. Indels Passing Filtering
        8. Sites Masked
        9. Coverage
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Minimum Coverage Allowed : ${minCov}" >> analysis-parameters.txt
    echo "Reference Supplied : ${refName}" >> analysis-parameters.txt
    echo "Primers Supplied : ${primerFileName}" >> analysis-parameters.txt
    echo "Medaka Model : ${model}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "Sample,Raw Reads,Mapped Reads,Clipped Mapped Reads,Average Read Depth,SNPs,Indels,Masked Sites,Coverage" > stats-summary.csv
    """
}

// Generates an alignment of the asembly contigs to a reference genome
// using minimap.
process MiniMap2_Alignment {
    input:
        // Tuple contains the sample basename
        // and the assembly directory
        tuple val(base), file(fastq)
        // The name of the output directory
        val outDir
        // Tuple contains the reference basename 
        // and the reference fasta file to be aligned to.
        tuple val(refName), file(ref) 

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file("${base}-align.bam")
        // A string containing statistics generated during this step
        env summary

    publishDir "${outDir}", mode: 'copy'
    
    script:
    /*
    Uses minimap2 to align reads to a reference.

    Then, samtools is used to convert the alignment sam into bam format
    (samtools view). The bam format is then sorted and stored in a bam
    file (samtools sort).

    The number of mapped reads is calculated using samtools and
    appended to the string containing previously collected statistics.
    */
    """
    #!/bin/bash

    raw_reads=\$((\$(cat ${fastq} | wc -l)/4)) 

    minimap2 -ax map-ont ${ref} ${fastq} > align.sam

    samtools view -b align.sam | samtools sort > ${base}-align.bam

    mapped_reads=\$(samtools view -F 0x04 -c ${base}-align.bam)

    summary="${base},\$raw_reads,\$mapped_reads"
    """
}

// Soft clips xGen primers from an alignment using primerclip
process XGen_Primer_Clip {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(bam)
        // The primer masterfile for input into primerclip.
        file primersFile
        // The name of the output directory
        val outDir
        // The number of threads provided
        val threads
        // A string containing existing statistics string 
        // to be added to.
        val existingSummary

    output:
        // Tuple contains the file basename and the clipped, sorted bam alignment file.
        tuple val(base), file("${base}-clipped-sorted.bam")
        // A string containing statistics generated during this step along with previously generated statistics
        env summary

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Primerclip takes a name sorted alignment file in sam format. Thus, the first
    step uses samtools to sort and convert the alignment bam.

    Then, primerclip is run to soft clip primers from the alignment.

    Finally, the primerclip output is converted back to bam format and sorted using
    samtools.

    The number of mapped reads after clipping is calculated and added 
    to the string containing previously collected statistics.
    */
    """
    #!/bin/bash

    samtools view -bS ${bam} | samtools sort -@ ${threads} -n -O sam > ${base}-sort.sam

    primerclip -s ${primersFile} ${base}-sort.sam ${base}-clip.sam

    samtools view -b ${base}-clip.sam | samtools sort -@ ${threads} > ${base}-clipped-sorted.bam

    clipped_mapped_reads=\$(samtools view -F 0x04 -c ${base}-clipped-sorted.bam)

    average_read_depth=\$(samtools depth -a -J -q 0 -Q 0 ${base}-clipped-sorted.bam | awk -F'\t' 'BEGIN{totalCov=0} {totalCov+=\$3} END{print totalCov/NR}')

    summary="${existingSummary},\$clipped_mapped_reads,\$average_read_depth"
    """
}

// Uses the 'consensus' step of the Medaka polishing tool
// to correct the alignment and produce an hdf file that
// can be used for subsequenct variant calling.
process Medaka_Consensus {
    input:
        // Tuple contains the sample base name and the clipped bam file.
        tuple val(base), file(bam)
        // The medaka model to use (supplied by the user)
        val model
        // The output directory
        val outDir
        // A string containing previously generated statistics to be added to.
        val existingSummary

    output:
        // Tuple contains the sample name, the corrected consensus by medaka,
        // and the clipped bam file.
        tuple val(base), file("${base}.hdf"), file(bam)
        // A string containing statistics generated during this step 
        // along with previously generated statistics
        env summary

    publishDir "${outDir}", mode:'copy'

    script:
    /*
    Uses samtools to index the bam file

    Then, runs the cosnensus step of the medaka polisher to create a
    polished hdf file from a bam.
    */
    """
    #!/bin/bash

    samtools index ${bam}

    medaka consensus ${bam} --model ${model} ${base}.hdf --debug

    summary="${existingSummary}"
    """
}

// Calls and filters variants from the polished alignment file
process Call_Variants {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(hdf), file(bam)
        // The script base directory name (to call python scripts)
        val baseDir
        // The output directory name
        val outDir
        // Tuple contains the reference name and reference fasta file
        tuple val(refName), file(ref)
        // The minimum coverage cutoff
        val minCov
        // A string containing previously generated statistics to be appended
        // to.
        val existingSummary

    output:
        // Tuple contains the file basename, alignment bamfile, filtered snp vcf, and filtered indel vcf
        tuple val(base), file(bam), file("${base}-snps-filtered.vcf"), file("${base}-indels-filtered.vcf")
        // Tuple contains the unfiltered vcf and multiallelic filtered vcf files.
        file "${base}-longshot.vcf"
        // A string containing statistics generated during this step 
        // along with previously generated statistics
        env summary

    publishDir "${outDir}", mode: 'copy'

    // NOTE: Have yet to encounter multiallelic site - once one is encountered, I can create
    // a script to convert these sites into biallelic sites.

    /*
    Indexes the bam file using samtools and then calls variants
    using medaka which takes both the alignment bam and the previously
    generated hdf file to avoid errors when making calls. Then, the vcf file
    is indexed using tabix and gzipped.
    
    However, a caveat is that medaka does not give enough statistcs in its 
    VCF file to filter baed on depth and alternate allele abundance. Thus,
    we use the tool Longshot, which will take a set of variants in VCF format
    and output them with more detail. We supply the following options to remove any default 
    filtering options from Longshot:

        - -P 0 - removes a p value cutoff for displaying variants based on
                 strand bias
        
        - -A - allows a larger maximum coverage when considering variants
        
        - -a 0 - sets the minimum base quality to 0 (allows for all bases to be considered
                  when calculating allele fractions)

        - --min_mapq 0 - sets the minimm mapping quality to 0 (allows all reads in the alignment
                         to be considered when calculating depth/allele fraction)
        
        - --no_haps - Longshot phases variants by default, this turns that option off

    
    Next, VCFTools is used to separate indels and SNPs into
    separate files which are indexed and gzipped.

    BCFTools is then used to filter the variants in each
    file. Variants at sites that fall below the minimum coverage are removed.
    As well, variants that are not the majority 
    (variant_allele_count/total_bases_at_site > reference_allele_count/total_bases_at_site).
    
    The number of passing SNPs and Indels are calculated and appended to the
    summary data.
    */
    script:
    """
    #!/bin/bash

    samtools index ${bam}

    medaka variant ${ref} ${hdf} ${base}-medaka.vcf

    bgzip ${base}-medaka.vcf
    tabix ${base}-medaka.vcf.gz

    longshot -P 0 -A -a 0 --min_mapq 0 --no_haps -v ${base}-medaka.vcf.gz --bam ${bam} --ref ${ref} --out ${base}-longshot.vcf

    vcftools --keep-only-indels --vcf ${base}-longshot.vcf --recode --recode-INFO-all --stdout > ${base}-indels.vcf
    bgzip ${base}-indels.vcf
    tabix ${base}-indels.vcf.gz

    vcftools --remove-indels --vcf ${base}-longshot.vcf --recode --recode-INFO-all --stdout > ${base}-snps.vcf
    bgzip ${base}-snps.vcf
    tabix ${base}-snps.vcf.gz

    bcftools view -i "((INFO/AC[0] + INFO/AC[1]) >= ${minCov}) && ((INFO/AC[1] / INFO/DP) > (INFO/AC[0] / INFO/DP))" ${base}-indels.vcf.gz > ${base}-indels-filtered.vcf
    num_indels=\$(grep -v "^#" ${base}-indels-filtered.vcf | wc -l)

    bcftools view -i "((INFO/AC[0] + INFO/AC[1]) >= ${minCov}) && ((INFO/AC[1] / INFO/DP) > (INFO/AC[0] / INFO/DP))" ${base}-snps.vcf.gz > ${base}-snps-filtered.vcf
    num_snps=\$(grep -v "^#" ${base}-snps-filtered.vcf | wc -l)

    summary="${existingSummary},\$num_snps,\$num_indels"
    """
}

// Generates a consensus genome by applying variants to the
// reference and masking based on read depth
process Generate_Consensus_Genome {
    input:
        // Tuple contains the file basename, the alignment bam, the snp vcf file
        // and the indel vcf file
        tuple val(base), file(bam), file(snps), file(indels)
        // The name of the base directory
        val baseDir
        // The name of the base directory
        val outDir
        // Tuple contains the reference file name and reference file
        tuple val(refName), file(ref)
        // The minimum coverage threshold
        val minCov
        // A string containing previously generated statistics to be appended to.
        val existingSummary

    output:
        // Tuple contains the consensus fasta and the sites that were masked in 
        // a bed file.
        tuple val(base), file("${base}-consensus.fasta")
        // The bed file containing the masked sites.
        file "${base}-mask-sites.bed"
        // A string containing statistics generated during this step 
        // along with previously generated statistics
        env summary

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    The script first computes sites to mask by identifying sites that have less than
    the minimum coverage provided. However, there are issues with the pileup
    format that make this difficult. Samtools mpileup's output has the depth in 0-based
    format, which makes it impossible to distinguish between a site with 0 and 1 coverage.
    
    Thus, the pipeline instead makes use of bedtools subtract. It first creates a pileup for only sites with
    coverage, uses an in-house script to filter sites with less than the minimum coverage,
    and converts these into a bed file.

    Next, the pipeline creates a pileup containing every site regardless of whether there is coverage
    and converts this into a bed file using the same in-house script. 

    Finally, the sites we want to keep (those above the minimum coverage threshold) are substracted
    from the bed file with every site, to give us the low-coverage sites.

    An interesting case occurs when no reads aligned to the provided reference. In this case, the pileup
    will be blank (regardless if it is told to output all sites), and if the process is followed, will
    give us no coordinates to mask (ending in 100% coverage). This is accounted for by checking if the
    bed file created from the pileup containing all sites is blank. If it is not blank, the pipeline continues normally.
    If it is, then a bed file containing every coordinate in the genome is created (as we would want the whole genome
    masked in this case).

    Additionally, there is another interesting case when the sites that fall within a deletion are marked as masked.
    Because masking is applied first, this will cause an error when applying the variants (as the deletion site will contain 
    an N character and will not match the VCF reference). Thus, the pipeline uses bcftools to create a bed file for all indel sites,
    and then subtracts these from the low coverage sites. Now, these results are the final sites to be masked.

    The bedtools maskfasta command is then used to mask the reference at these positions.

    Then the variants are applied to the mask fasta. The reason this is done after masking is 
    because the pileup (and therefore masking) positions do not account for indels, which would
    shift the genomic coordinates. If we then applied masking after applied indels, we would end 
    up masking regions that we should not have (the coordinates of these regions is now shifted).

    Finally, the fasta is wrapped to make it visually appealing. 

    The number of masked sites and genome coverage are calculated and appended to the 
    other statistics.
    */
    """
    #!/bin/bash

    samtools mpileup --no-BAQ -d 100000 -x -A -Q 0 -f ${ref} ${bam} > ${base}.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i ${base}.pileup -o passed-sites.bed --minCov ${minCov}

    samtools mpileup --no-BAQ -d 100000 -x -A -a -Q 0 -f ${ref} ${bam} > all-sites.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i all-sites.pileup -o all-sites.bed 

    if [[ -s all-sites.bed ]]; then
        bedtools subtract -a all-sites.bed -b passed-sites.bed > ${base}-low-cov-sites.bed

        bcftools query -f'%CHROM\t%POS0\t%END\n' ${indels} > indel-sites.bed

        bedtools subtract -a ${base}-low-cov-sites.bed -b indel-sites.bed > ${base}-mask-sites.bed
    else
        bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${ref} > ${base}-mask-sites.bed
    fi

    num_mask=\$(bioawk -c bed 'BEGIN{SITES=0} {SITES+=\$end-\$start } END{print SITES}' ${base}-mask-sites.bed)

    bedtools maskfasta -fi ${ref} -bed ${base}-mask-sites.bed -fo masked.fasta

    bgzip ${snps}
    tabix ${snps}.gz

    bgzip ${indels}
    tabix ${indels}.gz

    bcftools consensus -f masked.fasta ${snps}.gz > with-snps.fasta

    bcftools consensus -f with-snps.fasta ${indels}.gz > with-indels-snps.fasta

    bioawk -c fastx '{ gsub(/\\n/,"",seq); print ">${base}"; print \$seq }' with-indels-snps.fasta > ${base}-consensus.fasta

    seq_len=\$(bioawk -c fastx '{ print length(\$seq) }' < ${base}-consensus.fasta)

    coverage=\$(awk -v mask=\$num_mask -v len=\$seq_len 'BEGIN { print (1 - (mask / len)) * 100 }')

    summary="${existingSummary},\$num_mask,\$coverage"
    """
}


// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // Tuple contains the sample basename and forward/reverse reads (the basename
        // is the only value important to this function).
        val summary
        // The output directory.
        val outDir

    script:
    /*
    The summary string containing the statistics collected as the pipeline
    was run are appended to the summary file.
    */
    """
    #!/bin/bash

    echo "${summary}" >> ${outDir}/stats-summary.csv
    """  

}