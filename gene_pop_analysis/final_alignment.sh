#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=mafft_mp                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=7gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=kdr_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
pwd; hostname; date
export OMP_NUM_THREADS=1
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
module load slurm/18.08.4 shared DefaultModules Workspace/v1
module load ks575/Mafft/v7.471
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
#
base=$WORKSPACE/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis
work_dir=$base/all_consensus_fa/mod_headers
while read -r f
do
	cat ${work_dir}/${f}/*.fa >> ${work_dir}/${f}/${f}.all.fasta
	mafft --auto ${work_dir}/${f}/${f}.all.fasta > ${work_dir}/${f}/${f}.all.align.fasta
	cp -vrf ${work_dir}/${f}/${f}.all.fasta ${base}/final_fasta_files
	cp -vrf ${work_dir}/${f}/${f}.all.align.fasta ${base}/final_aligned_files
done < "${work_dir}/all.folders.name"
