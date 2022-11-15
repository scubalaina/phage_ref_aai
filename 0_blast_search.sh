#!/bin/sh

#SBATCH -J bee_inphared_blastn
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alainarw@vt.edu
#SBATCH --nodes=1
#SBATCH -t 100:00:00

##SBATCH -p normal_q

#SBATCH -A aylwardlab
#SBATCH --partition=aylward_lab

cd /groups/Aylward_Lab/Alaina/Side_hustles/bee_phage/beephage_inphared/


blastn -query rd2_checkvbee_orfs.fna -out rd2_checkvbee_inphared_blast.txt -db /groups/Aylward_Lab/riley/caudo_inphared/caudo_inphared_genes_blastdb -evalue 0.001 -outfmt 6