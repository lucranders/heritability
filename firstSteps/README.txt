The '.r' files presented here are called by the '.sh' files with the same name.
Every program has a number as a suffix, naming the order in which the files were processed.

The pipeline is described as following:
        
        * data_prep_1.sh: For each dataset ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz (in which ${chr} represents a number, ranging from 1 to 22, relative to each chromossome):

        	* Only samples which the hla expressions are known are filtered (bcftools) and
        	* Only alleles with at least one copy are kept (bcftools) then
        	* The list of alleles with correlation lesser than 0.1^(0.5) is generated (calling data_prep_1.r)

        