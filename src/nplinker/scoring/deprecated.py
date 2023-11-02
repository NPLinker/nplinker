def expand_spectrum_score(spectrum, gcf, scoring_function, strain_list):
    initial_score, initial_metadata = scoring_function(spectrum, gcf, strain_list)
    expanded_score, expanded_metadata = scoring_function(spectrum.family, gcf, strain_list)
    print(
        "{} <-> {}\tInitial: {}, expanded: {} ({} spectra in family (id = {}))".format(
            spectrum,
            gcf,
            initial_score,
            expanded_score,
            len(spectrum.family.spectra),
            spectrum.family.family_id,
        )
    )


# TODO needs updating due to annotation changes
def knownclusterblast_scoring(spectral_like, gcf_like, mibig_map):
    score = 0
    metadata = None
    if spectral_like.empty_default_annotations():
        print("No annotations")
        return None, None
    kcb = []
    for bgc in gcf_like.bgc_list:
        these = bgc.known_cluster_blast
        # if hasattr(bgc,'metadata'):
        #     these = bgc.metadata.get('knownclusterblast',None)
        if these:
            for mibig, score in these:
                kcb.append((mibig, score))
    if len(kcb) == 0:
        return None, None
    total_score = 0
    for annotation in spectral_like.get_annotations():
        for mibig, score in kcb:
            short_mibig = mibig.split("_")[0]
            if short_mibig in mibig_map:
                m = match(annotation, mibig_map[short_mibig])
                if m:
                    metadata = m
                    total_score += int(score)
                    print(m)
    return total_score, metadata
