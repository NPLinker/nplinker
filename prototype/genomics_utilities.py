from Bio import SeqIO
def get_known_cluster_blast(bgc):
	if hasattr(bgc,'antismash_file'):		
		bgc_file = bgc.antismash_file
		if not bgc_file:
			print("No antismash file in {}".format(bgc))
			return [] # Stops it being repeatedly called
		records = list(SeqIO.parse(bgc.antismash_file,'genbank'))
		hits = []
		for record in records:
			for feature in record.features:
				if feature.type == 'cluster':
					bgc_list = feature.qualifiers.get('knownclusterblast',None)
					if bgc_list:
						for bgc_item in bgc_list:
							tokens = bgc_item.split('\t')
							bgc_name = tokens[0].split('.')[1].strip()
							similarity = tokens[1].split('(')[1].split('%')[0]
							hits.append((bgc_name,similarity))
		return hits
