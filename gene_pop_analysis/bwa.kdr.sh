#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=HiC_Mp5                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=30                           # Number of CPU cores per task
#SBATCH --mem=5gb                                       # Job memory request
#SBATCH --time=250:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=p4505_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
pwd; hostname; date
export OMP_NUM_THREADS=30
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
module load shared DefaultModules gcc/8.2.0 slurm/18.08.4 Workspace/v1 
module load bwa/0.7.17 
module load samtools/1.9
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
#
cd /nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/USA
out_dir="/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/USA"
work_dir="/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/USA"
while IFS=$'\t' read -r f1 f2
do
	outfile=$(echo $f1 | cut -d'/' -f11 | cut -d'-' -f1)
	ref_dir="/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/indexes2"
	OIFS=$IFS; IFS=$'\n'; LINES=($(<$ref_dir/all.fasta.txt)); IFS=$OIFS
	for LINE in "${LINES[@]}"
	do
		name=$(echo $LINE | cut -d'/' -f12 | cut -d'.' -f2 | cut -d'_' -f2)
		echo "Started the sample $f1 and reference $name"
		bwa mem -t 30 ${ref_dir}/${name} "$f1" "$f2" | samtools view -bS -F 4 - | samtools sort -@ 30 -m 1G -o $out_dir/${outfile}.${name}.bam
	done
done < "${work_dir}/US.file.txt"
date
