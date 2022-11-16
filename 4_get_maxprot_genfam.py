import os, sys, re, statistics, argparse
from collections import defaultdict

args_parser = argparse.ArgumentParser(description="Script for calculating statistics on families each phage matched.", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-i', '--infile', required=True, help='Input table of ANI or AAI from step 2, with each phage and its AAI or ANI of each reference the phage hit to at least one protein.')
args_parser.add_argument('-o', '--outfile', required=True, help='Output file with family statistics.')

args_parser = args_parser.parse_args()

infile_a = args_parser.infile
outfile_a = args_parser.outfile

blast_ani = open(infile_a,'r')
outfile_open = open(outfile_a,'a')

genhit2ident = defaultdict(list)
genhit2prop = defaultdict(list)
genhit2prot = defaultdict(list)

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
		outfo = genome + "\t" + max_gen + "\t" + hinfo
		outfile_open.write(outfo + "\n")
