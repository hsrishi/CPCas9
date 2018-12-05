Order of scripts
	1. trim_reads_forBWA.py : trims/QCs NGS reads for use with BWA
	2. bwa_samtools.py : aligns trimmed NGS reads from (1) to premade dCas9 index and yields alignment mapping files
	3. process_SAM.py : processes alignment output from (2) and produces csv w/ nucleotide mapping of CP variants
	4. prep_data_for_DESeq.ipynb : converts output from (3) into a csv file w/ amino acid mapping of CP variants
	5. calculate_enrichment_DESeq.R : calculates enrichment scores for CP variants using mapping from (4)
	