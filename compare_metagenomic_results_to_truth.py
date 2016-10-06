"""
Compares count-based output of metagenomic analysis tools with known ground truth of dataset
Lorian Schaeffer, Pachter Lab, UC Berkeley
"""

import sys
import csv
import numpy
import string
import math
import itertools
import collections
import cPickle
import os
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'
import urllib
import urllib2
import argparse
import copy
import pprint
import scipy.stats

scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))
numpy.set_printoptions(precision=4)

class TaxaDict():
	# Consists of a dict mapping names -> taxids, and a dict mapping taxids -> Taxonomy objects
	def __init__(self, names = dict(), taxids = dict()):
		self.names = names
		self.taxids = taxids


	def get_taxa(self,name_or_id):
		taxa = self.get_taxa_by_id(name_or_id)
		if not taxa:
			taxa = self.get_taxa_by_name(name_or_id)
		return taxa

	def get_taxa_by_name(self,name):
		# get Taxonomy object by name
		try:
			return self.taxids[self.names[name]]
		except:
			try:
				return self.taxids[self.names[name.lower()]]
			except:
				return ''

	def get_taxa_by_id(self,taxid):
		# get Taxonomy object by ID
		try:
			return self.taxids[taxid]
		except:
			return ''

	def get_name_by_id(self,taxid):
		try:
			return self.taxids[taxid].name
		except:
			return ''

	def add_taxa(self,name,official_name,taxid,lineage):
		tax_entry = Taxonomy(official_name,taxid,lineage)
		self.names[str(name)] = taxid
		self.taxids[taxid] = tax_entry
		return tax_entry

	def add_blank(self):
		tax_entry = Taxonomy('Bacteria',2,{'superkingdom': '2', 'no rank': '131567'})
		self.names[str('bacteria')] = 2
		self.taxids[taxid] = tax_entry
		return tax_entry

	def remove_taxa(self,taxid):
		del self.names[self.get_taxa_by_id(taxid).name]
		del self.taxids[taxid]

# global taxa dictionary
tax_dict = TaxaDict()

class Taxonomy():

	def __init__(self, name = "", taxid = 0, lineage = dict()):
		self.name = name
		self.taxid = taxid
		self.lineage = lineage

	def __str__(self):
		return "{} | {} | {}".format(self.name,self.taxid,self.lineage)

	def filter_tier(self,rank):
		# return only if at or below a given taxonomic rank
		ranks = ['strain','species','genus','family','order','class','phylum']
		valid_ranks = set(ranks[0:ranks.index(rank)+1])
		if valid_ranks & set(self.lineage.keys()):
			return self.taxid
		else:
			return 0

	def filter_tier_exact(self,rank):
		# return only if AT a given taxonomic rank
		if rank in self.lineage.keys():
			return self.taxid
		else:
			return 0

	def return_tier(self,rank):
		# return the taxid for the specified rank
		try:
			return self.lineage[rank]
		except:
			return self.taxid

	def has_lineage(self,taxid):
		# does this Taxonomy object have this taxid somewhere in its lineage?
		return taxid in self.lineage.values()

class Dataset():

	def __init__(self, species=[], abundance=[], counts=[], size=[]):
		self.species = list(species)
		self.abundance = list(abundance)
		self.counts = list(counts)
		self.size = list(size) # can be empty
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts)) # defaults to 0

	def add_record(self, new_species, new_abundance, new_counts, new_size=None):
		self.species.extend([new_species])
		self.abundance.append(new_abundance)
		self.counts.append(new_counts)
		self.size.append(new_size)
		self.counts_by_sp[new_species] = new_counts

	def set_by_array(self, array):
		# Warning: this overwrites all existing data in this dataset
		try:
			self.species, self.abundance, self.counts, self.size = zip(*array)
		except:
			self.species, self.abundance, self.counts = zip(*array)
		self.species = list(self.species)
		self.abundance = [float(x) for x in self.abundance]
		self.counts = [float(x) for x in self.counts] # some programs like eXpress return counts as floats
		self.size = [int(x) for x in self.size]
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts))

	def remake_index(self):
		# when the species listing has updated and the dictionary needs to match
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts))

	def print_record(self, species):
		print "{0}:\n\t{1} ({2} reads)".format(species,self.lookup_abundance(species),self.lookup_count(species))

	def get_array(self):
		return zip(self.species,self.abundance,self.counts)

	def lookup_abundance(self, species):
		try:
			return self.abundance[self.species.index(species)]
		except:
			for sp in species.split('?'):
				try:
					return self.abundance[self.species.index(sp)]
				except:
					pass
		return 0

	def lookup_count(self, species):
		count = self.counts_by_sp[species]
		if count == 0:
			try:
				count = self.counts_by_sp[tax_dict.get_taxa_by_name(species).taxid]
			except:
				pass
		return count

	def lookup_size(self, species):
		try:
			return self.size[self.species.index(species)]
		except:
			return 0

	def null_list(self):
		# For use when one of the variables needs to be filled in
		return [0]*len(self.species)

	def clean_names(self):
		# Force clean_names to be lowercase and replace problematic symbols
		self.species = [s.lower().translate(string.maketrans("()[]:-","      ")).strip().replace(" ","_").replace(".","") for s in self.species]

	def sort_by_name(self):
		# Alphabetical by species
		self.set_by_array(sorted(self.get_array(),key=lambda x:x[0]))

	def set_threshold(self, threshold=1):
		# Drop entries with very low estimated counts
		self.set_by_array([r for r in self.get_array() if r[2] >= threshold])
		print "Number of non-zero abundances: {0}".format(len(self.species))

	def convert_to_percentage(self):
		# Normalize TPM/FPKM abundances
		total_ab = math.fsum(self.abundance)
		if total_ab == 0: # give up if there are no abundances
			return

		unknown = self.lookup_abundance('unknown')
		self.abundance = [100*v/(total_ab - unknown) for v in self.abundance]
		if unknown != 0: # prevent unknowns from being normalized
			self.abundance[self.species.index("unknown")] = unknown

		assert round(math.fsum(self.abundance)) == 100 + round(unknown)

	def remove_matches(self, target):
		self.set_by_array([r for r in self.get_array() if (r[0].find(target) == -1)])
		print "Number of entries after removing {0}: {1}".format(target,len(self.species))

	def delete_record(self,target):
		try:
			index = self.species.index(target)
		except:
			print "{} not found to delete".format(target)
			return
		del self.species[index]
		del self.counts[index]
		del self.abundance[index]
		del self.counts_by_sp[index]
		try:
			del self.size[index]
		except:
			return


def print_names(strains):
	print [tax_dict.get_name_by_id(s) for s in strains]

def uptier_taxa(strains,rank):
	uptier_names = filter_by_taxa(strains,rank) # gives entries of rank or below

	uptier_ids = []
	for s in strains:
		if s in uptier_names:
			rank_id = tax_dict.get_taxa_by_id(s).return_tier(rank)
			uptier_ids.append(rank_id)
		else:
			uptier_ids.append(s)

	lookup_tax_list(uptier_ids)
	return uptier_ids

def collapse_strains(strains,rank):
	""" Group strains together by taxa """
	just_rank = uptier_taxa(strains.species,rank)

	rank_combo = [x for x in zip(just_rank,strains.abundance,strains.counts,strains.size)]
	rank_dict = {}
	for s,a,c,z in rank_combo:
		rank_dict.setdefault(s,[0,0,0]) # stores ab,counts,size in that order
		rank_dict[s][0] = rank_dict[s][0] + a
		rank_dict[s][1] = rank_dict[s][1] + c
		rank_dict[s][2] = rank_dict[s][2] + z

	a,c,z = zip(*rank_dict.values())
	j_rank = Dataset(rank_dict.keys(),a,c,z)

	return j_rank

def collapse_duplicates(raw_data):
	# Create dictionary of lists of duplicates
	dup_data = raw_data.get_array()
	set_sp = {}
	set_ab = {}
	set_co = {}
	set_sz = {}
	set_plasmids = {}
	for sp,ab,co in dup_data:
		if 'taxid' in sp: # retain useful information
			name = sp.rpartition('|')[0] # last segment is usually the original chromosome etc name
		else:
			name = sp.partition('_gi|')[0].partition('|')[0].partition('_gca')[0] #the prepended strain name
		set_sp.setdefault(name,[]).append(sp)
		set_ab.setdefault(name,[]).append(ab)
		set_co.setdefault(name,[]).append(co)
		set_sz.setdefault(name,[]).append(co)

	assert(set_ab.keys() == set_co.keys() == set_sp.keys())

	# New, clean dataset for data without duplicates
	undupe = Dataset()

	# Note: we include plasmids in the count total solely because i100 was simulated to include 1x plasmid coverage.
	for k,v in set_sp.items():
		if len(v) == 1: # just add record directly if it has no duplicates
			undupe.add_record(k,set_ab[k][0],set_co[k][0],set_sz[k][0])
		else: # sum counts and average abundances
			undupe.add_record(k,math.fsum(set_ab[k])/len(v),math.fsum(set_co[k]),math.fsum(set_sz[k]))

	print "Number of entries after combining duplicates: {0}".format(len(undupe.species))

	return undupe

def get_taxid(original_name):

	# problematic names
	known_names = dict([('bacterium ellin514 strain ellin514','Pedosphaera parvula Ellin514'),
						('bacteroidetes sp. f0058','Bacteroidetes oral taxon 274 str. F0058'),
						('baumannia cicadellinicola str. hc','Baumannia cicadellinicola str. Hc (Homalodisca coagulata)'),
						('bifidobacterium longum infantis atcc 55813','Bifidobacterium longum subsp. longum ATCC 55813'),
						('borrelia afzelii pko clone','Borrelia afzelii PKo'),
						('brucella abortus bv. 1 str. 9-941 chromosome','Brucella abortus bv. 1 str. 9-941'),
						('brucella abortus bv. 3 str.','Brucella abortus bv. 3'),
						('buchnera aphidicola str. bp','Buchnera aphidicola str. Bp (Baizongia pistaciae)'),
						('buchnera aphidicola str. sg','Buchnera aphidicola str. Sg (Schizaphis graminum)'),
						('buchnera aphidicola str. cc','Buchnera aphidicola BCc'),
						('buchnera aphidicola str. lsr1','Buchnera aphidicola str. LSR1 (Acyrthosiphon pisum)'),
						('candidatus hamiltonella defensa 5at','Candidatus Hamiltonella defensa 5AT (Acyrthosiphon pisum)'),
						('candidatus pelagibacter sp. htcc7211 1105874033148','Candidatus Pelagibacter sp. HTCC7211'),
						('candidatus ruthia magnifica str. cm','Candidatus Ruthia magnifica str. Cm (Calyptogena magnifica)'),
						('candidatus sulcia muelleri str. hc','Candidatus_Sulcia_muelleri_str._Hc_(Homalodisca_coagulata)'),
						('capnocytophaga sputigena atcc 33612 strain capno','Capnocytophaga sputigena ATCC 33612'),
						('clostridiales bacterium 1 7 47faa strain','Clostridiales bacterium 1_7_47FAA'),
						('clostridiales genomosp. bvab3 upii9-5','Mageeibacillus indolicus UPII9-5'),
						('clostridium difficile complete genome strain cf5','clostridium difficile cf5'),
						('clostridium difficile complete genome strain m120','clostridium difficile m120'),
						('enterococcus faecalis tx 0109','enterococcus faecalis tx0109'),
						('enterococcus faecalis tx 0411','enterococcus faecalis tx0411'),
						('enterococcus faecalis tx 0855','enterococcus faecalis tx0855'),
						('enterococcus faecalis tx 0860','enterococcus faecalis tx0860'),
						('enterococcus faecalis tx 2134','enterococcus faecalis tx2134'),
						('enterococcus faecalis tx 4248','enterococcus faecalis tx4248'),
						('escherichia coli o150:h5 se15','Escherichia coli SE15'),
						('lactobacillus brevis gravesensis atcc 27305','Lactobacillus brevis subsp. gravesensis ATCC 27305'),
						('lactobacillus reuteri sd2112 atcc 55730','Lactobacillus reuteri SD2112'),
						('lactobacillus rhamnosus gg atcc 53103','Lactobacillus rhamnosus GG'),
						('methanocaldococcus infernus me c','Methanocaldococcus infernus ME'),
						('nostoc azollae','\'Nostoc azollae\' 0708'),
						('providencia alcalifaciens ban1 integrating conjugative element icepalban1','Providencia alcalifaciens Ban1'),
						('salmonella enterica subsp. arizonae serovar 62:z4z23:--','Salmonella enterica subsp. arizonae serovar 62:z4,z23:-'),
						('salmonella enterica subsp. enterica serovar 4','Salmonella enterica subsp. enterica serovar 4,[5],12:i:- str. CVM23701'),
						('selenomonas sp. 67h29bp f0410','Selenomonas sp. oral taxon 149 str. 67H29BP'),
						('staphylococcus aureus aureus atcc baa-39','Staphylococcus aureus subsp. aureus ATCC BAA-39'),
						('staphylococcus aureus aureus mn8','Staphylococcus aureus subsp. aureus MN8'),
						('staphylococcus aureus aureus tch130/st-72','Staphylococcus aureus subsp. aureus TCH130'),
						('staphylococcus aureus aureus tch60','Staphylococcus aureus subsp. aureus TCH60'),
						('staphylococcus aureus aureus tch70','Staphylococcus aureus subsp. aureus TCH70'),
						('staphylococcus aureus aureus usa300 tch959','Staphylococcus aureus subsp. aureus USA300_TCH959'),
						('streptococcus gallolyticus tx20005','Streptococcus gallolyticus subsp. gallolyticus TX20005'),
						('streptococcus sp. 73h25ap f0408','Streptococcus sp. oral taxon 071 str. 73H25AP'),
						('streptomyces coelicolor a3','Streptomyces coelicolor A3(2)'),
						('wolbachia endosymbiont of drosophila willistoni tsc','Wolbachia endosymbiont of Drosophila willistoni TSC#14030-0811.24'),
						('xanthomonas campestris pv. vasculorum ncppb702','Xanthomonas vasicola pv. vasculorum NCPPB 702'),
						('bordetella bronchiseptica strain rb50','Bordetella bronchiseptica RB50'),
						('chlamydia trachomatis l2buch-1proctitis','Chlamydia trachomatis L2b/UCH-1/proctitis'),
						('ignicoccus hospitalis kin4i','Ignicoccus hospitalis strain KIN4/I'),
						('lawsonia intracellularis phemn1-00','Lawsonia intracellularis PHE/MN1-00'),
						('mycobacterium bovis af212297','Mycobacterium bovis AF2122/97'),
						('escherichia coli str. k-12 substr. mg1655star','Escherichia coli str. K-12 substr. MG1655'),
						('clostridium botulinum e1 str. \'bont e','Clostridium botulinum E1 str. \'BoNT E Beluga\''),
						('bacteria',2), # taxa with non-bacterial entries in NCBI database
						('bacillus',1386),
						('thermosipho',2420),
						('yersinia',629),
						('erwinia tasmaniensis et199','Erwinia tasmaniensis ET1/99'),
						('erwinia tasmaniensis et1 99','Erwinia tasmaniensis ET1/99'),
						('erwinia tasmaniensis et1/99','Erwinia tasmaniensis ET1/99'),
						("erwinia tasmaniensis strain et1/99",'Erwinia tasmaniensis ET1/99'),
						("erwinia tasmaniensis strain et199",'Erwinia tasmaniensis ET1/99'),
						('lawsonia intracellularis phemn1 00','Lawsonia intracellularis PHE/MN1-00'),
						('ignicoccus hospitalis kin4 i','Ignicoccus hospitalis KIN4/I'),
						])

	tax_entry = tax_dict.get_taxa_by_name(original_name) # skip the lookup if taxa was already found
	if tax_entry:
		return tax_entry

	if 'taxid' in original_name:
		taxid = original_name.partition('taxid|')[2].partition('|')[0]
		return taxid

	name = original_name.partition('|')[0].lower().replace('_',' ').strip()
	if name in known_names.keys(): # problematic names
		name = known_names[name]

	tax_entry = tax_dict.get_taxa_by_id(name) # skip the lookup if taxa was already found for alternate name
	if tax_entry:
		return tax_entry

	try:
		if type(name) == int or name.isdigit(): #if we've already passed a taxid instead of a text name
			return name

		# Clean name for URL use; replace marks that cause problems for NCBI lookup with spaces
		url_name = urllib.quote_plus(name.translate(string.maketrans("()[]:","     ").replace('_',' ').strip()))
		# Look up taxonomy ID
		handle = Entrez.esearch(db="taxonomy", term=url_name)
		records = Entrez.read(handle)
		taxid = records['IdList']
		if taxid:
			return taxid

		# lookup manually instead of through biopython
		page_taxa = urllib2.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={}'.format(url_name))
		for line in page_taxa:
			line = string.replace(line,"\t","")
			if line.startswith("<Id>"):
				taxid = line[4:-6]
				return taxid

		# look up through uniprot
		page_taxa = urllib2.urlopen('http://www.uniprot.org/taxonomy/?query={}'.format(url_name))
		for line in page_taxa:
			if line.startswith("</script>"):
				taxid = line.partition("tr id=\"")[2].partition("\"")[0]
				return taxid

		if 'gca' in original_name.lower() or '_str_' in original_name.lower():
			name = original_name.lower().partition('_gca')[0].partition('_str_')[0]
			url_name = urllib.quote_plus(name.translate(string.maketrans("()[]:","     ").replace('_',' ').strip()))
			handle = Entrez.esearch(db="taxonomy", term=url_name)
			records = Entrez.read(handle)
			taxid = records['IdList']
			return taxid

		print "\t Unable to find {} when looked up as {}".format(original_name,url_name)
		return ""

	except Exception as e: # because when a long string of name lookups errors in the middle, it hurts
		print e
		print original_name
		return ""


def lookup_tax(original_name):
	global tax_dict

	taxid = get_taxid(original_name)
	if not taxid or isinstance(taxid,Taxonomy): # if it's blank or is already a taxonomy
		return taxid

	try:
		# Get taxonomy for id
		handle = Entrez.efetch(db="taxonomy", id=taxid, mode="text", rettype="xml")
		taxon = Entrez.read(handle)[0] # only grab the first

		taxid = taxon["TaxId"]
		official_name = taxon['ScientificName']
		lineage = dict()
		try:
			for t in taxon["LineageEx"]:
				lineage[t['Rank']] = t['TaxId']
		except:
			pass # no parent taxonomy

		if taxon['Rank'] == 'no rank':
			if 'species' in lineage.keys():
				lineage['strain'] = taxid
			else:
				lineage[str(taxid)] = taxid
		else:
			lineage[taxon['Rank']] = taxid

		if official_name.startswith('[Eubacterium'): # special case handling
			lineage['genus'] = '1730'

		tax_entry = tax_dict.add_taxa(official_name,official_name,taxid,lineage)
		tax_entry = tax_dict.add_taxa(original_name,official_name,taxid,lineage)
		print "{} -> {}".format(original_name,official_name)

		if '2759' in lineage.values(): # catch eukaryotes
			print "\t\t\tWarning: {} appears to be a eukaryote!".format(original_name)
			print tax_entry

	except Exception as e: # because when a long string of name lookups errors in the middle, it hurts
		print e
		print original_name
		return ""

	return tax_entry

def lookup_tax_list(species_list):
	species_tax = []
	for name in species_list:
		if not(len(name) == 8 and name[0:2].isalpha() and name[3:7].isdigit()):
			# skip unlookupable entries like CJG34856

			try:
				add_to_list = lookup_tax(name).taxid
			except:
				try:
					add_to_list = lookup_tax(name.rpartition('_')[0]).taxid
				except:
					print "Unable to find {} -- complete failure".format(name)
					add_to_list = 2 # bacteria as default, so length/order of list doesn't change
		species_tax.append(add_to_list)
	return species_tax

def fix_transcript_names(species):
	cleanup = string.maketrans('-()+/','_____')
	transcript_names = [(sp.partition("|")[0].translate(cleanup,'.[]') +"_"+
		"..".join([d for d in sp.split(':') if d.isdigit() and len(d)>1])) for sp in species]

	return transcript_names


def process_input(filename,program,fragmented=False):
	""" Pull species names, abundance, and counts out of input file """

	suffix = filename.rpartition('.')[2]
	input_file = open(filename,'r')

	raw_est = Dataset()
	if program == 'express':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[1]
		raw_est.clean_names()
		raw_est.counts = [float(i) for i in zip(*input_data)[6]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[10]]
		raw_est.size = raw_est.null_list()

	elif program == 'clark':
		input_csv = csv.reader(input_file, 'excel')
		input_data = [r for r in input_csv]
		input_data = input_data[1:-1] #remove header row and unknown line at end
		raw_est.species = [r[0]+"|taxid|"+r[1]+"|" for r in input_data] # combine name and taxid
		raw_est.clean_names()
		raw_est.counts = [float(i) for i in zip(*input_data)[3]]
		raw_est.abundance = [float(i) if i != '-' else 0.0 for i in zip(*input_data)[5]]

	elif program == 'kraken':
		if not filename.endswith('.report'):
			if os.path.exists(filename +'.p'): # full kraken output is slow to process, so look for saved processed version
				print "Loading kraken input from pickled file..."
				est = cPickle.load(open(filename +'.p','rb'))
				est.set_threshold()
				est.species = lookup_tax_list(est.species) # from this point on, species are taxids
				est.remake_index()
				est = collapse_duplicates(est)
				est.remake_index()
				return est
			input_csv = csv.reader(input_file, 'excel-tab')
			input_counts = collections.Counter()
			for r in input_csv:
				input_counts.update([r[1].rpartition(';')[2]])
			input_table = [(k.rpartition(';')[2],v) for k,v in input_counts.iteritems()]
			cPickle.dump(input_table,open(filename +'.p','wb')) # save processed form
			raw_est.species = zip(*input_table)[0]
			raw_est.counts = [float(i) for i in zip(*input_table)[1]]
			raw_est.abundance = [0]*len(raw_est.species)
		else: # nice shiny report format from kraken, much easier to process
			input_csv = csv.reader(input_file, 'excel-tab')
			input_data = [r for r in input_csv]
			input_data = input_data[1:] #remove unclassified row
			raw_est.species = ["{}|taxid|{}|".format(r[5].strip(),r[4]) for r in input_data] # combine name and taxid
			raw_est.counts = [float(i) for i in zip(*input_data)[2]]
			raw_est.abundance = [float(i) for i in zip(*input_data)[0]]
			raw_est.size = raw_est.null_list()

	elif program == 'kallisto':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[0]
		raw_est.counts = [float(i) for i in zip(*input_data)[3]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[4]]
		raw_est.size = raw_est.null_list()

	elif program == 'gasic':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[0]
		raw_est.counts = [float(i) for i in zip(*input_data)[2]]
		raw_est.abundance = [0]*len(raw_est.species)
		raw_est.size = raw_est.null_list()

	elif program == 'bracken':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = ["{}|taxid|{}|".format(r[0],r[1]) for r in input_data] # combine name and taxid
		raw_est.counts = [float(i) for i in zip(*input_data)[5]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[6]]
		raw_est.size = raw_est.null_list()

	else:
		print "File is not supported input type"

	print "Number of raw entries: {0}".format(len(raw_est.species))

	raw_est.clean_names()
	raw_est.remake_index()
	raw_est.remove_matches('rna') # remove specific genes
	raw_est.remove_matches('gene_')
	est = collapse_duplicates(raw_est)

	est.convert_to_percentage()
	est.sort_by_name()
	est.set_threshold()

	if program == 'kraken' and not filename.endswith('.report'):
		print 'Saving kraken output to file...'
		cPickle.dump(est,open(filename +'.p','wb'))

	est.species = lookup_tax_list(est.species) # from this point on, species are taxids
	est.remake_index()
	est = collapse_duplicates(est)
	est.remake_index()

	return est

def dataset_truth(dataset='i100'):
	# Truth is in format "species",abundance,counts,genome size where undefined counts are 0

	truth = Dataset()
	i100_csv = [r for r in csv.reader(open(os.path.join(scriptdir,'i100_truth.csv'),'r'), 'excel')]

	truth.set_by_array(list(i100_csv))

	#truth.set_threshold()
	truth.clean_names()
	truth.sort_by_name()
	truth.species = lookup_tax_list(truth.species) # from this point on, species are taxids
	truth.remake_index()
	return truth

def filter_by_taxa(taxids,rank):
	filtered_names = set()
	for taxid in taxids:
		id_taxa = tax_dict.get_taxa_by_id(taxid)
		filtered_names.add(id_taxa.filter_tier(rank))

	filtered_names.discard(0) # 0 is returned if a name doesn't pass the rank filter

	# Discard specific entries that don't get caught automatically
	if rank == 'strain':
		filtered_names.discard('83333') # Escherichia coli K-12

	return filtered_names

def calc_counts_error(truth,est):
	adjusted_counts = numpy.zeros(len(truth.counts))
	adjusted_true_counts = numpy.zeros(len(est.counts))

	# estimated counts for true species
	adjusted_counts = [est.lookup_count(sp) for sp in truth.species]

	print "Total number of samples: {}".format(len(est.counts))
	print "% counts mis-assigned to non-truth taxa: {:.2f}%".format(100*(sum(est.counts) - sum(adjusted_counts))/sum(est.counts))
	print "% counts not assigned at all: {:.2f}% ({:.0f} counts assigned)".format(100*(sum(truth.counts) - sum(est.counts))/sum(truth.counts),sum(est.counts))

	# normalize counts to only ones that mapped at all:
	normalization_factor = sum(truth.counts)/sum(est.counts) # total counts/counts assigned
	normalized_counts = numpy.array(adjusted_counts)*normalization_factor

	diff_n = 100*(numpy.array(normalized_counts) - numpy.array(truth.counts))/numpy.array(truth.counts)
	print "Average relative error of normalized (by a factor of {:.2f}) assigned counts: {:.2f}%".format(normalization_factor,numpy.mean(abs(diff_n)))
	diff_sq_n = [d*d/10000 for d in diff_n]
	print "Relative root mean squared error of normalized assigned counts: {:.2f}%".format(numpy.mean(diff_sq_n) ** (0.5) * 100)

	return diff_n,normalized_counts,normalization_factor

def climb_tree(taxid,dataset):
	# Given a taxid and a Dataset, calculate the total counts that would end up in that taxid
	count_sum = 0.0
	for sp in dataset.species:
		if tax_dict.get_taxa_by_id(sp).has_lineage(taxid):
			count_sum += dataset.lookup_count(sp)
	return count_sum

def climb_tree_verbose(taxid,dataset):
	# Given a taxid and a Dataset, calculate the total counts that would end up in that taxid
	count_sum = 0.0
	for sp in dataset.species:
		if tax_dict.get_taxa_by_id(sp).has_lineage(taxid):
			count_sum += dataset.lookup_count(sp)
			print tax_dict.get_taxa_by_id(sp)
	return count_sum


def graph_error(truth, est, adjusted_abundance, diff, expname, tier, norm_factor, show_graphs, save_graphs, errors, program='kallisto'):
	# These imports are here so script can run on server w/out graphics
	import plotfunctions
	import matplotlib
	import seaborn

	# format program name
	if program == 'clark':
		program = 'CLARK'
	if program == 'bracken':
		program = 'Bracken'

	# Drop any entry that is above the targeted taxa level
	# sometimes the strain-level filter filters out species-level genomes, so add truth back in
	filtered_names = list(filter_by_taxa(est.species,tier).union(set(truth.species)))

	filtered_est = Dataset()
	filtered_est.set_by_array([(sp,errors[sp],est.lookup_count(sp)) for sp in filtered_names])
	# abundance field contains stdev errors if they exist

	truth.species = [tax_dict.get_name_by_id(s) for s in truth.species]
	truth.remake_index()
	filtered_est.species = [tax_dict.get_name_by_id(s) for s in filtered_est.species]
	filtered_est.counts = list(numpy.array(filtered_est.counts)*norm_factor)
	filtered_est.remake_index()

	all_species = filtered_est.species
	all_est = [filtered_est.counts[filtered_est.species.index(sp)] if sp in filtered_est.species else 0 for sp in filtered_est.species]
	all_true = [truth.counts[truth.species.index(sp)] if sp in truth.species else 0 for sp in filtered_est.species]

	all_diff = []
	for i,a in enumerate(all_est):
		try:
			all_diff.append(100*(a - all_true[i]) / max(a,all_true[i]))
		except ZeroDivisionError:
			all_diff.append(0) # if both the estimate and the actual are 0, we're good here

	# graph true abundances
	if len(filtered_est.species) == len(truth.species):
		xmax = len(truth.species)
		x = numpy.array(range(0,xmax))

		ab_filter = zip(truth.species,truth.counts,adjusted_abundance)
		ab_filter.sort( key=lambda x: x[1],reverse=True )
		ab_species,ab_true,ab_adjusted = zip(*ab_filter)

		plotfunctions.plot_setup_pre("{} estimated counts at {}-level"
			.format(program,tier), xlabels = ab_species,
			xticks = range(0,xmax), xrotation = -90, yaxislabel = 'Counts')

		plotfunctions.plot(x, ab_true, color='blue', label='True')
		plotfunctions.plot(x,ab_adjusted, color='red', label='Estimated', plot_type = 'scatter')
		matplotlib.pyplot.gca().set_ylim(bottom=0.) # set x axis at y=0
		if save_graphs:
			plotfunctions.plot_setup_post(save_file = expname +'_'+ tier +'_counts.png', show=show_graphs)
		else:
			plotfunctions.plot_setup_post(legend=False)

	# graph abundances for all species, not just the true ones, if above mean
	else:
		if len(est.species) - filtered_est.counts.count(0) > len(truth.species)*1.25:
				mean_ab = min(truth.counts)/10 #10% of lowest actual count in truth
				print "Filtering out any estimated results under {} counts".format(mean_ab)
		else:
			mean_ab = 1

		true = Dataset()
		true.species = all_species
		true.counts = all_true
		est = Dataset()
		est.species = all_species
		est.counts = all_est

		present = zip(all_species,all_true,all_est)
		present.sort(key=lambda x: x[1], reverse=True)

		present_species, present_true, present_est = zip(*present)
		present_errors = [0]*500

		present_true = []
		present_est = []
		present_species = []
		present_errors = []

		for i,sp in enumerate(filtered_est.species):
			if sp in truth.species or filtered_est.lookup_count(sp)>mean_ab: #only include species if it has a non-zero estimate or non-zero actual abundance
				present_species.append(sp.replace('_',' ')) # presentation names
				present_errors.append(filtered_est.lookup_abundance(sp))
				present_est.append(filtered_est.lookup_count(sp))
				present_true.append(truth.lookup_count(sp))

		xmax = len(present_species)
		x = numpy.array(range(0,xmax))

		# sort based on true counts, then sort false positives by est counts
		all_sort = zip(present_species,present_true,present_est,present_errors)
		all_sort.sort( key=lambda x: x[2],reverse=True )
		all_sort.sort( key=lambda x: x[1],reverse=True )
		#all_sort.sort(key=lambda x: x[0]) #sort by name for troubleshooting

		present_species,present_true,present_est,present_errors = zip(*all_sort)

		plotfunctions.plot_setup_pre(
			"{} estimated counts at {}-level"
			.format(program,tier), xticks = range(0,xmax),
			xrotation = -90, yaxislabel = 'Counts', xlabels = present_species)

		plotfunctions.plot(x, present_true, color='blue', label="True")
		plotfunctions.plot(x, present_est, plot_type = 'error', color='red', label="Estimated", fmt='o', yerr=present_errors)
		matplotlib.pyplot.gca().set_ylim(bottom=0.)
		if save_graphs:
			plotfunctions.plot_setup_post(save_file = expname +'_'+ tier +'_sigcounts.png', show=show_graphs)
		else:
			plotfunctions.plot_setup_post(legend=False)

	return

def graph_est(est, expname, tier, show_graphs, save_graphs, errors):
	# These imports are here so script can run on server w/out graphics
	import plotfunctions
	import matplotlib
	import seaborn

	# Drop any entry that is above the targeted taxa level
	filtered_names = list(filter_by_taxa(est.species,tier))

	if len(est.species) > 100:
		min_abundances = numpy.mean(est.counts)
		min_abundances = 1000
		print "Filtering out any estimated results under {} counts".format(min_abundances)
	else:
		min_abundances = 10
	filtered_est = Dataset()
	for sp in filtered_names:
		if est.lookup_count(sp) > min_abundances: # filter out results with less than min abundances
			filtered_est.add_record(sp,est.lookup_count(sp),est.lookup_count(sp),est.lookup_size(sp))

	present_species = [tax_dict.get_name_by_id(s) for s in filtered_est.species]
	present_est = list(filtered_est.abundance)

	# graph est abundances
	xmax = len(present_species)
	present_sp = [x.replace('_',' ') for x in present_species]
	x = numpy.array(range(0,xmax))

	# sort based on est abundances
	all_filter = zip(present_sp,present_est)
	all_filter.sort( key=lambda x: x[1],reverse=True )
	fil_sp,fil_est = zip(*all_filter)
	#fil_est = numpy.log(fil_est)

	print "Results that passed filter:"
	# calculate % for each result
	fil_per = [round(100*n/math.fsum(fil_est),1) for n in fil_est]
	all_display = zip(fil_sp,fil_est,fil_per)
	pprint.pprint(all_display)

	plotfunctions.plot_setup_pre(
		"{}-level estimated counts"
		.format(tier), xlabels = fil_sp, xticks = range(0,xmax),
		xrotation = -90, yaxislabel = 'Counts')

	plotfunctions.plot(x, fil_est, color='red', plot_type = 'scatter')
	matplotlib.pyplot.gca().set_ylim(bottom=0.)
	if save_graphs:
		plotfunctions.plot_setup_post(save_file = expname +'_'+ tier +'_estabundances.png', show=show_graphs)
	else:
		plotfunctions.plot_setup_post(legend=False)

	return


def main(argv=sys.argv):
	"""
	Command line usage: python compare_metagenomic_results_to_truth.py
						[filename] [program (default kallisto)] <-g show graphs?>
						<-s save graphs?> <--taxa level (defaults to all)>
	"""

	parser = argparse.ArgumentParser(description='Compare output of metagenomic analysis tools with ground truth of dataset')
	parser.add_argument('filename', help='Output file of metagenomic analysis tool')
	parser.add_argument('program', nargs='?', default='kallisto', help='Source program that created output. Valid options are: kallisto, kraken, clark, gasic, express. Defaults to kallisto.')

	parser.add_argument('-g', '--show-graphs', action='store_true', help='Display graphs of calculated errors')
	parser.add_argument('-s','--save-graphs', action='store_true', help='Save graphs of calculated errors to file')
	parser.add_argument('--taxa', default='all', help='Desired taxa level of analysis. Accepts one of: strain, species, genus, phylum. Defaults to all levels.')
	parser.add_argument('--dataset', default='i100', help='Dataset truth to be compared to. Accepts: i100, no_truth. Defaults to i100.')
	parser.add_argument('--bootstraps', default='', help='Directory containing .tsv files for kallisto bootstraps, to be converted into errors.')

	args = parser.parse_args()

	filename = args.filename
	exp_name = filename.rpartition('.')[0] # will be used for graph-naming purposes
	program = args.program
	show_graphs = args.show_graphs
	save_graphs = args.save_graphs
	dataset = args.dataset

	print "Running comparison on {}\n".format(filename)

	global tax_dict
	pickle_len = 0
	if os.path.exists(os.path.join(scriptdir,'species_taxonomy.pickle')): # stored taxonomy information from previous runs
		print "Loading taxonomy dict..."
		tax_dict = cPickle.load(open(os.path.join(scriptdir,'species_taxonomy.pickle'),'rb'))
		pickle_len = len(tax_dict.names.keys())

	truth = dataset_truth(dataset)
	estimated = process_input(filename,program)

	est_j_species = collapse_strains(estimated,'species')
	est_j_genus = collapse_strains(estimated,'genus')
	est_j_phylum = collapse_strains(estimated,'phylum')
	true_j_species = collapse_strains(truth,'species')
	true_j_genus = collapse_strains(truth,'genus')
	true_j_phylum = collapse_strains(truth,'phylum')

	bootstrap_counts = collections.defaultdict(int)
	if args.bootstraps:
		# process each file with process_input, then make summary stats
		bootstraps = [f for f in os.listdir(args.bootstraps) if f.endswith('.tsv')]
		bootstrap_ests = []
		for f in bootstraps:
			bootstrap_ests.append(process_input(os.path.join(args.bootstraps,f),program,truth))
		for s in estimated.species: # collect counts for each species from each bootstrap
			#bootstrap_counts[s] = numpy.std([b.lookup_count(s) for b in bootstrap_ests])
			bs = [b.lookup_count(s) for b in bootstrap_ests]
			ci = scipy.stats.t.interval(0.99, len(bs)-1, loc=numpy.mean(bs), scale=scipy.stats.sem(bs))
			bootstrap_counts[s] = int((ci[1] - ci[0])/2)

	if pickle_len != len(tax_dict.names.keys()): # don't re-pickle if nothing new is added
		print "Saving taxonomy dict..."
		cPickle.dump(tax_dict,open(os.path.join(scriptdir,'species_taxonomy.pickle'),'wb'))

	if program == 'clark' or program == 'bracken': # clark and bracken don't do strain-level assignment
		dataset_pairs = [('species',true_j_species,est_j_species),('genus',true_j_genus,est_j_genus),('phylum',true_j_phylum,est_j_phylum)]
	else:
		dataset_pairs = [('strain',truth,estimated),('species',true_j_species,est_j_species),('genus',true_j_genus,est_j_genus),('phylum',true_j_phylum,est_j_phylum)]

	if args.taxa == 'strain':
		dataset_pairs = [('strain',truth,estimated)]
	elif args.taxa == 'species':
		dataset_pairs = [('species',true_j_species,est_j_species)]
	elif args.taxa == 'genus':
		dataset_pairs = [('genus',true_j_genus,est_j_genus)]
	elif args.taxa == 'phylum':
		dataset_pairs = [('phylum',true_j_phylum,est_j_phylum)]

	for label,true,est in dataset_pairs:
		if dataset == 'no_truth':
			graph_est(est, exp_name, label, show_graphs, save_graphs, bootstrap_counts)

		else:
			print "\n{}-LEVEL ERROR:".format(label.upper())
			diff, adjusted_abundance, norm_factor = calc_counts_error(true,est)
			if show_graphs or save_graphs:
				graph_error(true, est, adjusted_abundance, diff, exp_name, label, norm_factor, show_graphs, save_graphs, bootstrap_counts, args.program)

if __name__ == "__main__":
	main()
