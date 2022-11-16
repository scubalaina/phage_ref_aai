import os, sys, re, shlex, subprocess, statistics, argparse
from collections import defaultdict
from itertools import combinations

args_parser = argparse.ArgumentParser(description="Script for calculating ani or aai from blast against a reference", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-b', '--blastout', required=True, help='Output from blast search')
args_parser.add_argument('-g', '--query_genes', required=True, help='Multi-fasta file of protein (amino acid) or gene/orfs (nucleotide) sequences of your phages.')
args_parser.add_argument('-r','--ref_genes', required=False, default = str(0.75), help='Multi-fasta file of protein (amino acid) or gene/orfs (nucleotide) sequences of REFERENCE phages.')
args_parser.add_argument('-t','--taxa_tbl', required=True, help='Tab-separated file with reference accessions in column 1 and family in column 3')
args_parser.add_argument('-o','--outfile', required=True, help='Name of output ANI/AAI table file')

args_parser = args_parser.parse_args()

blastout_a = args_parser.blastout
phage_seqs_a = args_parser.query_genes
ref_seqs_a = args_parser.ref_genes
ref_taxa_a = args_parser.taxa_tbl
outfile_a = args_parser.outfile



blast_open = open(sys.argv[1],'r')
phage_seqs = open(gene_seqs_a,'r')
ref_seqs = open(ref_seqs_a,'r')
ref_taxa = open(ref_taxa_a,'r')
outfile_open = open(outfile_a,'w')


gen2prots = defaultdict(int)
genome_list = []

for i in phage_seqs:
  line = i.rstrip()
  if line.startswith(">"):
    header = re.sub(">","",line)
    headerlist = header.split(" # ")
    prot = headerlist[0]
    protlist = prot.split("_")
    genlist = protlist[0:-1]
    genome = "_".join(genlist)
    #print(genome)
    gen2prots[genome] += 1
    genome_list.append(genome)

for i in ref_seqs:
  line = i.rstrip()
  if line.startswith(">"):
    header = re.sub(">","",line)
    headerlist = header.split(" # ")
    prot = headerlist[0]
    protlist = prot.split("_")
    genlist = protlist[0:-1]
    genome = "_".join(genlist)
    #print(genome)
    gen2prots[genome] += 1
    genome_list.append(genome)

gen2fam = {}

for i in ref_taxa:
  line = i.rstrip()
  tabs = line.split("\t")
  if line.startswith("Genome"):
    pass
  else:
    genome = tabs[0]
    family = tabs[2]
    gen2fam[genome] = family

protgen2ident = defaultdict(float)
prot2gen_info = defaultdict(lambda:"NA")
gen2hitgen = defaultdict(list)

for i in blast_open:
  line = i.rstrip()
  tabs = line.split("\t")
  prot = tabs[0]
  hit = tabs[1]
  ident = float(tabs[2])
  protlist = prot.split("_")
  genlist = protlist[0:-1]
  genome = "_".join(genlist)
  hitlist = hit.split("_")
  hitgenlist = hitlist[0:-1]
  hitgen = "_".join(hitgenlist)
  protgen = prot + "\t" + hitgen
  gen2hitgen[genome].append(hitgen)
  #print(hitgen)
  if genome == hitgen:
    pass
  else:
    if protgen in protgen2ident.keys():
      if ident > protgen2ident[protgen]:
        protgen2ident[protgen] = ident
      else:
        pass
    else:
      protgen2ident[protgen] = ident
      
genhit_ident = defaultdict(list)
genhit_count = defaultdict(list)
# python 1_get_inphared_ani.py rd2_checkvbee_inphared_blast.txt rd2_checkvbee_orfs.fna /groups/Aylward_Lab/riley/caudo_inphared/caudo_inphared_genes.fna /groups/Aylward_Lab/inphared/inphared_1Sep2021/1Sep2021_inphared_family_tbl.txt rd2_checkvbee_inphared_ani.tsv


for key, values in protgen2ident.items():
  protgen = key
  identity = values
  protgenlist = protgen.split("\t")
  prot = protgenlist[0]
  hitgen = protgenlist[1]
  protlist = prot.split("_")
  genlist = protlist[0:-1]
  genome = "_".join(genlist)
  genhit = genome + "\t" + hitgen
  genhit_count[genhit].append(prot)
  genhit_ident[genhit].append(identity)

genome_set = set(genome_list)
tbl_header_list = ["G1","G2","Family","ANI", "G1_hit_prots","Prop_G1_prot","G1_prots","G2_prots"]
tbl_header = "\t".join(tbl_header_list)
outfile_open.write(tbl_header + "\n")

for key, values in gen2hitgen.items():
  genome = key
  hit_gen_set = set(values)
  for i in hit_gen_set:
    genhit = genome + "\t" + i
    gen1 = genome
    gen2 = i
    hit_count_set = set(genhit_count[genhit])
    hit_count = len(hit_count_set)
    hit_ident_list = genhit_ident[genhit]
    if len(hit_ident_list) > 0:
        avg_ident = statistics.mean(hit_ident_list)
        prop_gen1 = hit_count / gen2prots[gen1]
        family = gen2fam[gen2]
        infolist = [gen1, gen2, family, str(avg_ident),str(hit_count),str(prop_gen1),str(gen2prots[gen1]),str(gen2prots[gen2])]
        info = "\t".join(infolist)
        outfile_open.write(info + "\n")
