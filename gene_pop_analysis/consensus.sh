#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=CYP_cons                         # Job name
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
module load bcftools/1.9
module load ks575/Tabix/v0.2.6
module load VCFtools/0.1.16
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
#
base=$WORKSPACE/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/Adam_P450_seq_analysis
work_dir=$base/all_CYP450s
ref_dir=$base/new_names
while IFS=$'\t' read -r f1
do
	outfile=$(echo $f1 | cut -d'/' -f12 | cut -d'.' -f1)
	name=$(echo $f1 | cut -d'/' -f12 | cut -d'.' -f2)
	echo "Started the sample $f1 and reference $name"
	bcftools mpileup -Ou -f ${ref_dir}/${name}.fasta "$f1" | bcftools call -Ou -mv | bcftools norm -f ${ref_dir}/${name}.fasta -Oz -o ${work_dir}/${outfile}.${name}.vcf.gz
	vcftools --gzvcf ${work_dir}/${outfile}.${name}.vcf.gz --minQ 40 --recode --recode-INFO-all --out ${work_dir}/${outfile}.${name}.fil
	bcftools view ${work_dir}/${outfile}.${name}.fil.recode.vcf -Oz -o ${work_dir}/${outfile}.${name}.fil.recode.vcf.gz
	tabix ${work_dir}/${outfile}.${name}.fil.recode.vcf.gz
	bcftools consensus -f ${ref_dir}/${name}.fasta ${work_dir}/${outfile}.${name}.fil.recode.vcf.gz > ${work_dir}/${outfile}.${name}.fa
done < "${base}/dedup.bam.list"
