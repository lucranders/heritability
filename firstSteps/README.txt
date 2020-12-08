The '.r' files presented here are called by the '.sh' files with the same name.
Every program has a number as a suffix, naming the order in which the files were processed.

The pipeline is described as following:
        
        * data_prep_1.sh: 

        	* For each dataset ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz (in which ${chr} represents a number, ranging from 1 to 22, relative to each chromossome):

        		* Only samples which the hla expressions are known are filtered (bcftools) AND
        		* Only alleles with at least one copy are kept (bcftools) - final file named as chr${chr}.vcf.gz - THEN
        		* The list of alleles with correlation lesser than 0.1^(0.5) is generated, named as snpSet_${chr}.rds (calling data_prep_1.r)

		* data_prep_2.sh: 

			* All the 22 datasets chr${chr}.vcf.gz are merged (bcftools) - named as allChr.vcf.gz - THEN
			* All lists of pruned alleles (generated in the previous step) are gathered into one vector structure and saved as fullPrunedList.rds AND
			* A structure to remove all unnecessary files based on regex is created, but not executed 