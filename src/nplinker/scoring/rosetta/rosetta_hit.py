class RosettaHit:
    def __init__(self, spec, gnps_id, mibig_id, bgc, spec_match_score, bgc_match_score):
        self.spec = spec
        self.gnps_id = gnps_id
        self.mibig_id = mibig_id
        self.bgc = bgc
        self.spec_match_score = spec_match_score
        self.bgc_match_score = bgc_match_score

    def __str__(self):
        return "RosettaHit: {}<-->{} via ({} ({:.3f}), {} ({:.3f}))".format(
            self.spec.id,
            self.bgc.id,
            self.gnps_id,
            self.spec_match_score,
            self.mibig_id,
            self.bgc_match_score,
        )

    def __repr__(self):
        return str(self)
