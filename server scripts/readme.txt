gatk.SNP.batch and gatk.wt.batch are the main script files. the others are for utility purposes and you might not need them

gatk.wt.batch assumes a single wt file and has defined input/outputs

gatk.SNP.batch is set up for a slurm batch array job, so you might need to change the array numbers at the start of the file to match your number of files