#!/bin/sh

#SBATCH -J phage_inphared_blastp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alainarw@vt.edu
#SBATCH --nodes=1
#SBATCH -t 100:00:00

##SBATCH -p normal_q

#SBATCH -A aylwardlab
#SBATCH --partition=aylward_lab

cd /groups/Aylward_Lab/Alaina/Side_hustles/bee_phage/beephage_inphared/


blastp -query your_phage_prots.faa -out your_phage_inphared_blast_prot.txt -db /groups/Aylward_Lab/inphared_nov22/inphared_11Nov2022/caudoviricetes/blastdbs/inphared_caudo_drep90_nov22_prot_blastdb -evalue 0.001 -outfmt 6
