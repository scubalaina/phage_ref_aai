import os, sys, re, statistics
from collections import defaultdict

genhit2ident = defaultdict(list)
genhit2prop = defaultdict(list)
genhit2prot = defaultdict(list)

blast_ani = open(sys.argv[1],'r')

gen2prots = defaultdict(list)
genprot2hit = defaultdict(list)

genhit2info = {}
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
		gen2prots[genome].append(prots)
		genprot2hit[genome + "\t" + str(prots)].append(hit)
		hinfo = [fam,str(ani),str(prots),str(prop)]
		genhit = genome + "\t" + hit
		genhit2info[genhit] = hinfo

for key, values in gen2prots.items():
	genome = key
	maxprot = max(values)
	maxhits = genprot2hit[genome + "\t" + str(maxprot)]
	if len(maxhits) > 1:
		pass
	else:
		max_gen = "".join(maxhits)
		hinfolist = genhit2info[genome + "\t" + max_gen]
		hinfo = "\t".join(hinfolist)
		print(genome + "\t" + max_gen + "\t" + hinfo)

