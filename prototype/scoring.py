def metcalf_scoring(spectral_like,gcf_like,strains,both = 10,met_not_gcf = -10,gcf_not_met = 0,neither = 1):
	cum_score = 0
	for strain in strains:
		in_spec = spectral_like.has_strain(strain)
		in_gcf = gcf_like.has_strain(strain)
		if in_spec and in_gcf:
			cum_score += both
		if in_spec and not in_gcf:
			cum_score += met_not_gcf
		if in_gcf and not in_spec:
			cum_score += gcf_not_met
		if not in_gcf and not in_spec:
			cum_score += 1
	return cum_score