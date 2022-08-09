#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=QC                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=30                              # Number of CPU cores per task
#SBATCH --mem=1gb                                       # Job memory request
#SBATCH --time=600:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=QC_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=30
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load slurm/18.08.4 shared Workspace/v1 samtools/1.9 bwa/0.7.17
cd /nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/QC_bams
list=/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis
work_dir=/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis/QC_bams
OIFS=$IFS; IFS=$'\n'; LINES=($(<$list/bam.list)); IFS=$OIFS
for LINE in "${LINES[@]}"
do
	name=`echo $LINE | cut -d'/' -f12 | cut -d'.' -f1-2`
	samtools sort -@ 16 -m 1G -n -o $work_dir/$name.readname.bam $LINE
	samtools fixmate -m -@ 30 $work_dir/$name.readname.bam $work_dir/$name.fixmate.bam
	samtools sort -@ 30 -m 1G -o $work_dir/$name.fx.sorted.bam $work_dir/$name.fixmate.bam
	samtools markdup -s -@ 30 $work_dir/$name.fx.sorted.bam $work_dir/$name.dedup.bam
	samtools index -@ 30 $work_dir/$name.dedup.bam
	rm -vrf $work_dir/$name.readname.bam $work_dir/$name.fixmate.bam $work_dir/$name.fx.sorted.bam
done
