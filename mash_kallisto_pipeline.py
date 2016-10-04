'''
The mash results are tab delimited lists of Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes:
Bacillus_aryabhattai_gca_001043825.GCA_001043825.1.30.dna.genome.fa     /home/lorian/scratch/illumina_100species_trimmed.1.fq   1       1       0/1000
Bacillus_atrophaeus_1942.GCA_000165925.1.30.dna.genome.fa       /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.295981        0.000944763     1/1000
Bacillus_bombysepticus_str_wang.GCA_000831065.1.30.dna.genome.fa        /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.197292        2.3647e-28      8/1000
'''

# run from inside directory with ALL genomes; will move genomes one level up OR to specified --directory
import csv
import argparse
import collections
import os
import re
import urllib2
from Bio import Entrez
Entrez.email = ''

def count_sp(fastas):
	# pick out unique species, more or less
	unique_species = collections.Counter()
	for f in fastas:
		name = f.partition('GCA')[0].partition('gca')[0].partition('.')[0].strip('_').strip('.')
		name_parts = name.split('_')
		if len(name_parts) == 2: # just the species is left
			sp_name = name
		elif 'endosymbiont' in name_parts: # hard to parse these
			sp_name = name
		elif '_sp._' in name or '_sp_' in name or '_species_' in name:
			sp_name = '_'.join(name_parts[0:3])
		else:
			sp_name = '_'.join(name_parts[0:2])

		if not (f.strip('_').startswith(sp_name+"_") or f.strip('_').startswith(sp_name+".")): # catch weird shit
			print "{} -> {} -> {}??".format(f,name,sp_name)

		unique_species.update([sp_name])
	return unique_species

def get_taxid(uid):
	# look up taxid
	if uid.startswith('GCA_'):
		page_taxa = urllib2.urlopen('http://www.ebi.ac.uk/ena/data/view/{}&display=xml'.format(uid))
		for line in page_taxa:
			if line.strip().startswith('<TAXON_ID>'):
				return line.partition('>')[2].partition('<')[0]

	if uid.startswith('NC_') or uid.startswith('gi_'):
		if uid.startswith('gi_'):
			uid = uid.partition('gi_')[2] #NCBI wants them as just straight numbers
		handle = Entrez.efetch("nucleotide", id=uid, retmode="xml")
		records = Entrez.read(handle)
		# standard location
		taxid = records[0]['GBSeq_feature-table'][0]['GBFeature_quals'][-1]['GBQualifier_value'].partition(':')[2]
		if not taxid:
			# sometimes it's not the last quals list, so we have to search for it
			taxid = [r['GBQualifier_value'] for r in records[0]['GBSeq_feature-table'][0]['GBFeature_quals'] if ('taxon' in r['GBQualifier_value'])][0].partition(':')[2]
		return taxid

	return ''

def collapse_contigs(f):
	basename = f.partition('.dna.genome.fa')[0].partition('.mfa')[0]
	try:
		mfa = open(f,'r')
	except IOError:
		print "{} does not exist".format(f)
		return False

	text = ""

	firstline = True
	for line in mfa:
		if firstline:
			# replace header with name
			new_name = basename.partition('.1.30')[0].partition('.GCA')[0].replace(" ", "_")
			# try to filter out only the useful parts of the original header: GCA, NC, or gi uids
			name_parts = filter(None, re.split("[|:]+", line.replace('gi|','gi_').replace('_gi', '|gi')))
			keep_ids = [w for w in name_parts if (w.startswith('GCA_') or w.startswith('NC_') or w.startswith('gi_'))]
			for id in keep_ids:
				taxid = get_taxid(id)
				if taxid:
					break
			if not taxid:
				print "Unable to lookup taxid for {}".format(f)
				print line
				print keep_ids
			
			text = ">{0}|kraken:taxid|{1}|{2}\n".format(new_name,taxid,"|".join(keep_ids))
			firstline = False
		elif line.startswith('>'):
			text += 'NNNNNNNNNN' #indicate possible gaps between chr, plasmids, shotgun pieces, etc
		elif line!= "\n": #skip empty lines within fasta
			text += line

	text = text + "\n" #add newline to end of file for eventual cat
	fa = open(basename + '.cat.fa', 'w')
	fa.seek(0)
	fa.write(text)
	fa.truncate()
	fa.close()
	return basename + '.cat.fa'

parser = argparse.ArgumentParser(description='Process mash results for kallisto index creation. Run script from within directory containing all genome fasta files.')
parser.add_argument('filename', help='Mash output file')
parser.add_argument('top_strains', help="How many strains of each species to keep for the quantification step")
parser.add_argument('--directory', default="../", help="Directory to put files for kallisto index creation. Default is moving them one directory up.")
parser.add_argument('--dry-run', action="store_true", help="If set, lists files that would be moved, but does not create or move any files.")
args = parser.parse_args()

with open(args.filename,'r') as mash_file:
	mash_csv = csv.reader(mash_file, delimiter='\t')
	mash_data = [r for r in mash_csv]
	mash_1 = dict()
	mash_5 = dict()
	truth = dict()
	for r in mash_data:
		matching_hashes = int(r[4].split('/')[0])
		if r[0][0].islower():
			truth[r[0]] = matching_hashes
		if matching_hashes > 0:
			mash_1[r[0]] = matching_hashes
		if matching_hashes > 4:
			mash_5[r[0]] = matching_hashes

species = count_sp(mash_5.keys())
sp_map = {sp.lower(): [(st,mash_1[st]) for st in mash_5.keys() if sp.lower() in st.lower()] for sp in species}
final_st = []

for sp in sp_map.keys():
	# pick top N of each species
	sp_map[sp].sort(key=lambda x:x[1],reverse=True)
	final_st.extend(sp_map[sp][:min(int(args.top_strains),len(sp_map[sp]))])

final_names = list(set(zip(*final_st)[0]))
final_names.sort()
with open('mash_names.txt','w') as f:
	f.write('\n'.join(final_names) + '\n')
print "Files to be moved listed in mash_names.txt."

if not args.dry_run: #don't create or move files
	for f in final_names:
		new_name = collapse_contigs(f)
		
		if new_name:
			os.system("mv {0} {1}".format(new_name,args.directory))
			if not os.path.exists(os.path.join('{0}','{1}').format(args.directory,new_name)):
				print new_name

	#compare lists
	moved_files = [f for f in os.listdir(args.directory) if f.endswith(".cat.fa")]
	for f in final_names:
		basename = f.partition('.dna.genome.fa')[0].partition('.mfa')[0]
		if not basename+'.cat.fa' in moved_files:
			print "Missing: {}".format(f)
		
	print "All hits: {}\t At least 5 hits: {}\t Species: {}\t Strains kept: {}\t Strains moved: {}\t".format(len(mash_1.keys()), len(mash_5.keys()), len(species), len(final_names), len(moved_files))
		
else:
	print "All hits: {}\t At least 5 hits: {}\t Species: {}\t Strains kept: {}\t".format(len(mash_1.keys()), len(mash_5.keys()), len(species), len(final_names))
