#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name cont_scaffolds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o cont_scaffolds.%j.out
#SBATCH -A lp_svbelleghem

#Now it's necessary to study the synteny between the most contiguous (i.e. largest and longest) scaffolds. This is done by observing the AGP files from RagTag and extracting the largest
#scaffolds that match with the Fst signals from the popgen analysis

#To extract the linkage groups from the primary assemblies
awk '/^>LG01_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG01_dud.fasta
awk '/^>LG02_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG02_dud.fasta
awk '/^>LG03_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG03_dud.fasta
awk '/^>LG04_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG04_dud.fasta
awk '/^>LG05_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG05_dud.fasta
awk '/^>LG06_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG06_dud.fasta
awk '/^>LG07_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG07_dud.fasta
awk '/^>LG08_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG08_dud.fasta
awk '/^>LG09_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG09_dud.fasta
awk '/^>LG10_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG10_dud.fasta

awk '/^>LG01_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG01_nieu.fasta
awk '/^>LG02_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG02_nieu.fasta
awk '/^>LG03_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG03_nieu.fasta
awk '/^>LG04_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG04_nieu.fasta
awk '/^>LG05_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG05_nieu.fasta
awk '/^>LG06_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG06_nieu.fasta
awk '/^>LG07_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG07_nieu.fasta
awk '/^>LG08_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG08_nieu.fasta
awk '/^>LG09_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG09_nieu.fasta
awk '/^>LG10_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG10_nieu.fasta

#The to extract specific scaffolds (given the positions found in the corresponding AGP files). It's necessary to add a new heading since awk only extracts the raw sequence

#SCAFFOLD 1
awk 'NR==2 {print substr($0, 34117388, 72878678-34117388+1)}' LG01_dud.fasta > ext_scaffold_1_prim_dud.fa
sed -i '1i >scaffold_1_RagTag' ext_scaffold_1_prim_dud.fa

awk 'NR==2 {print substr($0, 1, 58882969-1+1)}' LG06_nieu.fasta > ext_scaffold_1_prim_nieu.fa
sed -i '1i >scaffold_1_RagTag' ext_scaffold_1_prim_nieu.fa

ntSynt ext_scaffold_1_prim_dud.fa ext_scaffold_1_prim_nieu.fa -p ext_scaffold_1 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_1.synteny_blocks.tsv --fai ext_scaffold_1_prim_dud.fa.fai ext_scaffold_1_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_1.synteny_blocks.tsv --sort_order ext_scaffold_1_prim_dud.fa.fai ext_scaffold_1_prim_nieu.fa.fai --fais > ext_scaffold_1.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_1_prim_dud.fa.fai ext_scaffold_1_prim_nieu.fa.fai --prefix ext_scaffold_1 --blocks ext_scaffold_1.synteny_blocks.sorted.tsv --length 100 --colour ext_scaffold_1_prim_dud.fa

cp ext_scaffold_1.links.tsv $VSC_DATA
cp ext_scaffold_1.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_1.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_1.links.sorted.tsv && mv ext_scaffold_1.links.sorted.tsv ext_scaffold_1.links.tsv
cat ext_scaffold_1.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_1.sequence_lengths.sorted.tsv && mv ext_scaffold_1.sequence_lengths.sorted.tsv ext_scaffold_1.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_1.sequence_lengths.tsv -l ext_scaffold_1.links.tsv --scale 25000000 --p ext_scaffold_1

#SCAFFOLD 2
awk 'NR==2 {print substr($0, 9100059, 39608137-9100059+1)}' LG01_dud.fasta > ext_scaffold_2_prim_dud.fa
sed -i '1i >scaffold_2_RagTag' ext_scaffold_2_prim_dud.fa

awk 'NR==2 {print substr($0, 43076330, 84460820-43076330+1)}' LG03_nieu.fasta > ext_scaffold_2_prim_nieu.fa
sed -i '1i >scaffold_2_RagTag' ext_scaffold_2_prim_nieu.fa

ntSynt ext_scaffold_2_prim_dud.fa ext_scaffold_2_prim_nieu.fa -p ext_scaffold_2 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_2.synteny_blocks.tsv --fai ext_scaffold_2_prim_dud.fa.fai ext_scaffold_2_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_2.synteny_blocks.tsv --sort_order ext_scaffold_2_prim_dud.fa.fai ext_scaffold_2_prim_nieu.fa.fai --fais > ext_scaffold_2.
synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_2_prim_dud.fa.fai ext_scaffold_2_prim_nieu.fa.fai --prefix ext_scaffold_2 --blocks ext_scaffold_2.synteny_blocks.sorted.tsv --l
ength 100 --colour ext_scaffold_2_prim_dud.fa

cp ext_scaffold_2.links.tsv $VSC_DATA
cp ext_scaffold_2.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_2.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_2.links.sorted.tsv && mv ext_scaffold_2.links.sorted.tsv ext_scaffold_2.links.tsv
cat ext_scaffold_2.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_2.sequence_lengths.sorted.tsv && mv ext_scaffold_2.sequence_lengths.sorted.tsv ext_scaffold_2.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_2.sequence_lengths.tsv -l ext_scaffold_2.links.tsv --scale 25000000 --p ext_scaffold_2

#SCAFFOLD 3
awk 'NR==2 {print substr($0, 3681807, 34037187-3681807+1)}' LG02_dud.fasta > ext_scaffold_3_prim_dud.fa
sed -i '1i >scaffold_3_RagTag' ext_scaffold_3_prim_dud.fa

awk 'NR==2 {print substr($0, 54672198, 93883070-54672198+1)}' LG01_nieu.fasta > ext_scaffold_3_prim_nieu.fa
sed -i '1i >scaffold_3_RagTag' ext_scaffold_3_prim_nieu.fa

ntSynt ext_scaffold_3_prim_dud.fa ext_scaffold_3_prim_nieu.fa -p ext_scaffold_3 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_3.synteny_blocks.tsv --fai ext_scaffold_3_prim_dud.fa.fai ext_scaffold_3_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_3.synteny_blocks.tsv --sort_order ext_scaffold_3_prim_dud.fa.fai ext_scaffold_3_prim_nieu.fa.fai --fais > ext_scaffold_3.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_3_prim_dud.fa.fai ext_scaffold_3_prim_nieu.fa.fai --prefix ext_scaffold_3 --blocks ext_scaffold_3.synteny_blocks.sorted.tsv --length 100 --colour ext_scaffold_3_prim_dud.fa

cp ext_scaffold_3.links.tsv $VSC_DATA
cp ext_scaffold_3.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_3.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_3.links.sorted.tsv && mv ext_scaffold_3.links.sorted.tsv ext_scaffold_3.links.tsv
cat ext_scaffold_3.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_3.sequence_lengths.sorted.tsv && mv ext_scaffold_3.sequence_lengths.sorted.tsv ext_scaffold_3.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_3.sequence_lengths.tsv -l ext_scaffold_3.links.tsv --scale 25000000 --p ext_scaffold_3

#SCAFFOLD 4
awk 'NR==2 {print substr($0, 74064473, 103390037-74064473+1)}' LG01_dud.fasta > ext_scaffold_4_prim_dud.fa
sed -i '1i >scaffold_4_RagTag' ext_scaffold_4_prim_dud.fa

awk 'NR==2 {print substr($0, 7103728, 45349787-7103728+1)}' LG02_nieu.fasta > ext_scaffold_4_prim_nieu.fa
sed -i '1i >scaffold_4_RagTag' ext_scaffold_4_prim_nieu.fa

ntSynt ext_scaffold_4_prim_dud.fa ext_scaffold_4_prim_nieu.fa -p ext_scaffold_4 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_4.synteny_blocks.tsv --fai ext_scaffold_4_prim_dud.fa.fai ext_scaffold_4_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_4.synteny_blocks.tsv --sort_order ext_scaffold_4_prim_dud.fa.fai ext_scaffold_4_prim_nieu.fa.fai --fais > ext_scaffold_4.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_4_prim_dud.fa.fai ext_scaffold_4_prim_nieu.fa.fai --prefix ext_scaffold_4 --blocks ext_scaffold_4.synteny_blocks.sorted.tsv --length 100 --colour ext_scaffold_4_prim_dud.fa

cp ext_scaffold_4.links.tsv $VSC_DATA
cp ext_scaffold_4.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_4.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_4.links.sorted.tsv && mv ext_scaffold_4.links.sorted.tsv ext_scaffold_4.links.tsv
cat ext_scaffold_4.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_4.sequence_lengths.sorted.tsv && mv ext_scaffold_4.sequence_lengths.sorted.tsv ext_scaffold_4.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_4.sequence_lengths.tsv -l ext_scaffold_4.links.tsv --scale 25000000 --p ext_scaffold_4

#SCAFFOLD 7
awk 'NR==2 {print substr($0, 22340581, 46000200-22340581+1)}' LG03_dud.fasta > ext_scaffold_7_prim_dud.fa
sed -i '1i >scaffold_7_RagTag' ext_scaffold_7_prim_dud.fa

awk 'NR==2 {print substr($0, 120572473, 153501276-120572473+1)}' LG01_nieu.fasta > ext_scaffold_7_prim_nieu.fa
sed -i '1i >scaffold_7_RagTag' ext_scaffold_7_prim_nieu.fa

ntSynt ext_scaffold_7_prim_dud.fa ext_scaffold_7_prim_nieu.fa -p ext_scaffold_7 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_7.synteny_blocks.tsv --fai ext_scaffold_7_prim_dud.fa.fai ext_scaffold_7_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_7.synteny_blocks.tsv --sort_order ext_scaffold_7_prim_dud.fa.fai ext_scaffold_7_prim_nieu.fa.fai --fais > ext_scaffold_7.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_7_prim_dud.fa.fai ext_scaffold_7_prim_nieu.fa.fai --prefix ext_scaffold_7 --blocks ext_scaffold_7.synteny_blocks.sorted.tsv --length 100 --colour ext_scaffold_7_prim_dud.fa

cp ext_scaffold_7.links.tsv $VSC_DATA
cp ext_scaffold_7.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_7.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_7.links.sorted.tsv && mv ext_scaffold_7.links.sorted.tsv ext_scaffold_7.links.tsv
cat ext_scaffold_7.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_7.sequence_lengths.sorted.tsv && mv ext_scaffold_7.sequence_lengths.sorted.tsv ext_scaffold_7.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_7.sequence_lengths.tsv -l ext_scaffold_7.links.tsv --scale 25000000 --p ext_scaffold_7

#SCAFFOLD 8
awk 'NR==2 {print substr($0, 1, 22185742-1+1)}' LG09_dud.fasta > ext_scaffold_8_prim_dud.fa
sed -i '1i >scaffold_8_RagTag' ext_scaffold_8_prim_dud.fa

awk 'NR==2 {print substr($0, 1, 32539133-1+1)}' LG09_nieu.fasta > ext_scaffold_8_prim_nieu.fa
sed -i '1i >scaffold_8_RagTag' ext_scaffold_8_prim_nieu.fa

ntSynt ext_scaffold_8_prim_dud.fa ext_scaffold_8_prim_nieu.fa -p ext_scaffold_8 -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv ext_scaffold_8.synteny_blocks.tsv --fai ext_scaffold_8_prim_dud.fa.fai ext_scaffold_8_prim_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks ext_scaffold_8.synteny_blocks.tsv --sort_order ext_scaffold_8_prim_dud.fa.fai ext_scaffold_8_prim_nieu.fa.fai --fais > ext_scaffold_8.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai ext_scaffold_8_prim_dud.fa.fai ext_scaffold_8_prim_nieu.fa.fai --prefix ext_scaffold_8 --blocks ext_scaffold_8.synteny_blocks.sorted.tsv --length 100 --colour ext_scaffold_8_prim_dud.fa

cp ext_scaffold_8.links.tsv $VSC_DATA
cp ext_scaffold_8.sequence_lengths.tsv $VSC_DATA

cat ext_scaffold_8.links.tsv  | mlr --tsv sort -f strand -n block_id > ext_scaffold_8.links.sorted.tsv && mv ext_scaffold_8.links.sorted.tsv ext_scaffold_8.links.tsv
cat ext_scaffold_8.sequence_lengths.tsv | mlr --tsv sort -f seq_id > ext_scaffold_8.sequence_lengths.sorted.tsv && mv ext_scaffold_8.sequence_lengths.sorted.tsv ext_scaffold_8.sequence_lengths.tsv
Rscript plot_synteny_blocks_gggenomes.R -s ext_scaffold_8.sequence_lengths.tsv -l ext_scaffold_8.links.tsv --scale 25000000 --p ext_scaffold_8

#SCAFFOLD 9

#SCAFFOLD 11 DUDZELE
awk 'NR==2 {print substr($0, 24993515, 44107433-24993515+1)}' LG05_dud.fasta > ext_scaffold_11_prim_dud.fa
sed -i '1i >scaffold_11_RagTag' ext_scaffold_11_prim_dud.fa

#SCAFFOLD 12 NIEUWPOORT
awk 'NR==2 {print substr($0, 24298258, 52146979-24298258+1)}' LG04_nieu.fasta > ext_scaffold_12_prim_nieu.fa
sed -i '1i >scaffold_12_RagTag' ext_scaffold_12_prim_nieu.fa

#Mapping for Asynt, map Dudzele and Nieuwpoort individual scaffolds against each other with minimap2 
conda activate thesis

minimap2 -x asm5 ext_scaffold_1_prim_dud.fa ext_scaffold_1_prim_nieu.fa | gzip > ext_scaffold_1.paf.gz
minimap2 -x asm5 ext_scaffold_2_prim_dud.fa ext_scaffold_2_prim_nieu.fa | gzip > ext_scaffold_2.paf.gz
minimap2 -x asm5 ext_scaffold_3_prim_dud.fa ext_scaffold_3_prim_nieu.fa | gzip > ext_scaffold_3.paf.gz
minimap2 -x asm5 ext_scaffold_4_prim_dud.fa ext_scaffold_4_prim_nieu.fa | gzip > ext_scaffold_4.paf.gz
minimap2 -x asm5 ext_scaffold_7_prim_dud.fa ext_scaffold_7_prim_nieu.fa | gzip > ext_scaffold_7.paf.gz
minimap2 -x asm5 ext_scaffold_8_prim_dud.fa ext_scaffold_8_prim_nieu.fa | gzip > ext_scaffold_8.paf.gz
#minimap2 -x asm5 ext_scaffold_9_prim_dud.fa ext_scaffold_9_prim_nieu.fa | gzip > ext_scaffold_9.paf.gz
#minimap2 -x asm5 ext_scaffold_11_prim_dud.fa ext_scaffold_11_prim_nieu.fa | gzip > ext_scaffold_11.paf.gz
