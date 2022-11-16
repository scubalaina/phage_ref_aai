import os, sys, re, statistics
from collections import defaultdict

gen_fam2genomes = defaultdict(int)
gen_fam2ident = defaultdict(list)
gen_fam2prop = defaultdict(list)
gen_fam2prot = defaultdict(list)

blast_ani = open(sys.argv[1],'r')

for i in blast_ani:
	line = i.rstrip()
	tabs = line.split("\t")
	if line.startswith("G1"):
		pass
	else:
		genome = tabs[0]
		hit = tabs[1]
		fam = tabs[2]
		ani = float(tabs[3])
		prots = int(tabs[4])
		prop = float(tabs[5])
		gen_fam = genome + "\t" + fam
		gen_fam2genomes[gen_fam] += 1
		gen_fam2ident[gen_fam].append(ani)
		gen_fam2prop[gen_fam].append(prop)
		gen_fam2prot[gen_fam].append(prots)

headerlist = ["Genome","Family","Number_genomes_hit", "Max_ani","Max_prots","Max_prop"]
print("\t".join(headerlist))

for key, values in gen_fam2genomes.items():
	keylist = key.split("\t")
	genome = keylist[0]
	fam = keylist[1]
	avg_ani = max(gen_fam2ident[key])
	avg_prop = max(gen_fam2prop[key])
	avg_prot = max(gen_fam2prot[key])
	num_gen = str(values)
	infolist = [key, num_gen, str(avg_ani), str(avg_prot), str(avg_prop)]
	print("\t".join(infolist))
