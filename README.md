# GenoCaller_indent
GenoCaller_indent is a python program that calls genotypes from bam files at positions/regions specified in a bed file while ignoring the some specified number of base pairs at either end of the read (i.e. the indent). This is primarily for use with UDG-half treated ancient DNA (Rohland et al. 2015)

This version of the program will take any genotype with a low quality heterozygote call (Q<30) and convert to the next best homozygote call. To turn this off, edit the python script by changing line 43 to "GQ=0". pysam and matplotlib libraries in your version of python need to be installed to run the script.

Three files are created: An emit all vcf file noting the call for all base pairs in the bed file. A vcf file with only those sites with evidence for at least one alternative allele and that passes a predetermined QUAL filter. A haploid emit all vcf (ish), that gives you the most likely base under a haploid model. If two or more basepairs are tied, the reported allele is randomly chosen.

To run aDNA_GenoCaller type:

./GenoCaller_indent.py indexed_bamfile bed_file reference_genome indent

An additional script is included that can be run in the same way (GenoCaller_indent_randread.py) which performs genotyping using only a single random read from a target site. The output is an emit all vcf.
