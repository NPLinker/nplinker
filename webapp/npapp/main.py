import os

from bokeh.models.widgets import RadioGroup, Slider, Div, Select, CheckboxGroup
from bokeh.models.widgets import Button, TextInput, PreText, RadioButtonGroup
from bokeh.models import CustomJS
from bokeh.plotting import curdoc
from bokeh.models import ColumnDataSource

from nplinker.metabolomics import Spectrum, MolecularFamily
from nplinker.genomics import GCF, BGC
from nplinker.annotations import gnps_url, GNPS_KEY
from nplinker.scoring.methods import LinkCollection

from searching import SEARCH_OPTIONS, Searcher

from tables_functions import NA_ID

from tooltips import create_popover, wrap_popover
from tooltips import TOOLTIP_TABLES_FILTERING, TOOLTIP_HYBRID_SCORING

GENOMICS_SCORING_MODES = ['manual', 'GCF']
METABOLOMICS_SCORING_MODES = ['manual', 'MolFam']

GENOMICS_MODE_BGC, GENOMICS_MODE_GCF = range(2)
METABOLOMICS_MODE_SPEC, METABOLOMICS_MODE_MOLFAM = range(2)

# Scoring modes are:
# 1. BGCs (i.e. the immediate parent GCFs) against Spectra
# 2. GCFs (i.e. all GCFs containing selected BGCs) against Spectra
# 3. BGCs (...) against MolFams
# 4. GCFs (...) against MolFams
# 5. Spectra against GCFs
# 6. MolFams against GCFs
# (there are 6 and not 8 because BGCs aren't "first-class" objects as far as scoring is concerned)

# for display/output purposes, 1+3 and 2+4 are basically identical, because you still end
# up having a list of GCFs as input to the scoring function

SCO_MODE_BGC_SPEC, SCO_MODE_BGC_MOLFAM,\
    SCO_MODE_GCF_SPEC, SCO_MODE_GCF_MOLFAM,\
    SCO_MODE_SPEC_GCF, SCO_MODE_MOLFAM_GCF = range(6)

SCO_MODE_NAMES = ['BGCs to Spectra', 'BGCs to MolFams', 'GCFs to Spectra', 'GCFs to MolFams',
                'Spectra to GCFs', 'MolFams to GCFs']

SCO_GEN_TO_MET = [True, True, True, True, False, False]

# basic template for a bootstrap "card" element to display info about a
# given GCF/Spectrum etc
# format params:
# - unique id for header, e.g spec_heading_1
# - colour for header
# - unique id for body, e.g. spec_body_1
# - title text
# - id for body, same as #2
# - parent id 
# - main text
TMPL = open(os.path.join(os.path.dirname(__file__), 'templates/tmpl.basic.py.html')).read()

# same as above but with an extra pair of parameters for the "onclick" handler and button id
TMPL_ON_CLICK = open(os.path.join(os.path.dirname(__file__), 'templates/tmpl.onclick.py.html')).read()

TMPL_SEARCH = open(os.path.join(os.path.dirname(__file__), 'templates/tmpl.basic.search.py.html')).read()

TMPL_SEARCH_ON_CLICK = open(os.path.join(os.path.dirname(__file__), 'templates/tmpl.onclick.search.py.html')).read()

# what gets displayed in the <div> for results when there are no results
RESULTS_PLACEHOLDER = """
    <div class="result-placeholder"><h3>No results to display</h3></div>
"""

class ScoringHelper(object):
    """
    This class is a wrapper around some of the details of handling the various
    different scoring modes supported by the webapp
    """

    GENOMICS_LOOKUP  = {GENOMICS_MODE_BGC: [SCO_MODE_BGC_SPEC, SCO_MODE_BGC_MOLFAM],
                        GENOMICS_MODE_GCF: [SCO_MODE_GCF_SPEC, SCO_MODE_GCF_MOLFAM]}

    METABOLOMICS_LOOKUP = [SCO_MODE_SPEC_GCF, SCO_MODE_MOLFAM_GCF]

    def __init__(self, nh):
        self.method_and_mode = True
        self.genomics_mode = GENOMICS_MODE_BGC
        self.metabolomics_mode = METABOLOMICS_MODE_SPEC
        self.mode = SCO_MODE_BGC_SPEC
        # True if scoring is currently being done *from* genomics side *to* metabolomics
        self.from_genomics = True
        self.nh = nh
        self.clear()

    def clear(self):
        self.bgcs = []
        self.gcfs = []
        self.spectra = []
        self.molfams = []
        self.inputs = set()
        self.scoring_results = LinkCollection()
    
    def set_genomics(self):
        self.from_genomics = True
        self._update()

    def set_metabolomics(self):
        self.from_genomics = False
        self._update()

    def update_genomics(self, m):
        self.from_genomics = True
        self.genomics_mode = m
        self._update()

    def update_metabolomics(self, m):
        self.from_genomics = False
        self.metabolomics_mode = m
        self._update()

    def gen_to_met(self):
        return SCO_GEN_TO_MET[self.mode]

    def _update(self):
        if self.from_genomics:
            self.mode = ScoringHelper.GENOMICS_LOOKUP[self.genomics_mode][self.metabolomics_mode]
        else:
            self.mode = ScoringHelper.METABOLOMICS_LOOKUP[self.metabolomics_mode]

    def set_results(self, scoring_results):
        self.scoring_results = scoring_results

    @property
    def mode_name(self):
        return SCO_MODE_NAMES[self.mode]

    def generate_scoring_objects(self, bgcs, spectra, tsne_id):
        """
        This method handles constructing the set of objects to be used as input to
        the nplinker "get_links" scoring method, based on the type of objects selected
        from the plots and the currently selected TSNE projection (if any, currently not used).

        For most scoring modes it doesn't need to do anything much, but for the
        GCF to Spec and GCF to MolFam modes, the initial selection of BGCs 
        may reference GCFs that don't appear in the current projection 
        (as in these modes we are saying "select ALL GCFs that contain ANY of 
        the selected BGCs"). In this case the unavailable GCFs need to be filtered
        out. 
        """
        scoring_objs = set()

        # store these for later use
        self.bgcs = bgcs
        self.spectra = spectra

        if self.mode == SCO_MODE_BGC_SPEC or self.mode == SCO_MODE_BGC_MOLFAM:
            # in either of these two modes, we simply want to select the set 
            # of "parent" GCFs for each of the selected BGCs, which is simple
            # (although have to filter out dups here!)
            scoring_objs = set([parent for bgc in bgcs for parent in bgc.parents])
            self.gcfs = list(scoring_objs)
        elif self.mode == SCO_MODE_GCF_SPEC or self.mode == SCO_MODE_GCF_MOLFAM:
            # this pair of modes are more complex. first build a set of *all* 
            # GCFs that contain *any* of the selected BGCs

            # this is the set of GCFs available in the current TSNE projection
            available_gcfs = self.nh.available_gcfs[tsne_id]
            all_bgcs = set()

            # for each selected BGC...
            for bgc in bgcs:
                if bgc not in self.nh.bgc_gcf_lookup:
                    print('WARNING: missing BGC in lookup: {}'.format(bgc))
                    continue

                # ... get the set of GCFs of which it is a member...
                gcfs_for_bgc = self.nh.bgc_gcf_lookup[bgc]
                print('BGC {} parents = {}'.format(bgc, gcfs_for_bgc))

                # ... and filter out any not in the current TSNE projection
                gcfs_for_bgc = [gcf for gcf in gcfs_for_bgc if gcf in available_gcfs]

                scoring_objs.update(gcfs_for_bgc)

                for gcf in gcfs_for_bgc:
                    all_bgcs.update(gcf.bgcs)

            print('Selection of {} BGCs => {} GCFs'.format(len(bgcs), len(scoring_objs)))
            self.gcfs = list(scoring_objs)
            # now need to update the set of BGCs selected to include all of those
            print('Final total of BGCs: {}'.format(len(all_bgcs)))
            self.bgcs = list(all_bgcs)
        elif self.mode == SCO_MODE_SPEC_GCF:
            # in this mode can simply use the list of spectra
            scoring_objs = spectra
        elif self.mode == SCO_MODE_MOLFAM_GCF:
            # TODO
            # here i think will need to go from Spectra -> MolFams in same way as
            # for BGCs -> GCFs above
            print('NOT DONE YET')
            return None
        else:
            raise Exception('Unknown scoring mode! {}'.format(self.mode))

        return list(scoring_objs)

class NPLinkerBokeh(object):

    def __init__(self, helper, table_session_data):
        # handle the "helper" object from the server_lifecycle module, which deals
        # with loading the dataset and creating data ready to be wrapped in DataSources
        self.nh = helper
        self.table_session_data = table_session_data

        self.score_helper = ScoringHelper(self.nh)

        # by default filter out scoring results where the pair of objects has no shared strains
        self.filter_no_shared_strains = True

        self.searcher = Searcher(self.nh.nplinker)

    @property
    def nplinker(self):
        return self.nh.nplinker

    def generate_gcf_spec_result(self, pgindex, gcf, link_data):
        # generate a GCF >> Spectra nested list

        body = '<div class="accordion" id="accordion_gcf_{}">'.format(pgindex)
        # generate the body content of the card by iterating over the spectra it matched
        j = 0
        
        # sort by metcalf scoring if possible 
        # TODO this will need some more thought when other methods added!
        if self.scoring_objects['metcalf'] in self.score_helper.scoring_results.methods:
            sorted_links = self.score_helper.scoring_results.get_sorted_links(self.scoring_objects['metcalf'], gcf)
        else:
            sorted_links = link_data.values()

        for link in sorted_links:
            spec = link.target
            score = '[' + ', '.join('{}={}'.format(m.name, m.format_data(link[m])) for m in link.methods) + ']'

            spec_hdr_id = 'spec_result_header_{}_{}'.format(pgindex, j)
            spec_body_id = 'spec_body_{}_{}'.format(pgindex, j)
            spec_title = 'Spectrum(parent_mass={:.4f}, id={}), scores <strong>{}</strong>, shared strains=<strong>{}</strong>'.format(spec.parent_mz, spec.spectrum_id, score, len(link.shared_strains))
            if spec.has_annotations():
                spec_title += ', # annotations={}'.format(len(spec.annotations))

            spec_body = self.generate_spec_info(spec, link)

            # set up the chemdoodle plot so it appears when the entry is expanded 
            spec_btn_id = 'spec_btn_{}_{}'.format(pgindex, j)
            spec_plot_id = 'spec_plot_{}_{}'.format(pgindex, j)
            spec_body += '<center><canvas id="{}"></canvas></center>'.format(spec_plot_id)
            span_badges = '&nbsp|&nbsp;'.join(['<span class="badge badge-info">{}</span>'.format(m.name) for m in link.methods])

            # note annoying escaping required here, TODO better way of doing this?
            spec_onclick = 'setupPlot(\'{}\', \'{}\', \'{}\');'.format(spec_btn_id, spec_plot_id, spec.to_jcamp_str())

            hdr_color = 'ffe0b5'
            if len(spec.annotations) > 0:
                hdr_color = 'ffb5e0'
            body += TMPL_ON_CLICK.format(hdr_id=spec_hdr_id, hdr_color=hdr_color, btn_target=spec_body_id, btn_onclick=spec_onclick, btn_id=spec_btn_id, 
                                         btn_text=spec_title, span_badges=span_badges, body_id=spec_body_id, body_parent='accordion_gcf_{}'.format(pgindex), body_body=spec_body)
            j += 1

        body += '</div>'
        return body

    def generate_spec_gcf_result(self, pgindex, spec, link_data):
        # generate a Spectrum >> GCF nested list

        body = '<div class="accordion" id="accordion_spec_{}>'.format(pgindex)
        # now generate the body content of the card by iterating over the GCFs it matched
        
        j = 0

        # sort by metcalf scoring if possible 
        # TODO this will need some more thought when other methods added!
        if self.scoring_objects['metcalf'] in self.score_helper.scoring_results.methods:
            sorted_links = self.score_helper.scoring_results.get_sorted_links(self.scoring_objects['metcalf'], spec)
        else:
            sorted_links = link_data.values()

        for link in sorted_links:
            gcf = link.target
            score = '[' + ', '.join('{}={}'.format(m.name, m.format_data(link[m])) for m in link.methods) + ']'

            # construct a similar object info card header for this GCF
            gcf_hdr_id = 'gcf_result_header_{}_{}'.format(pgindex, j)
            gcf_body_id = 'gcf_body_{}_{}'.format(pgindex, j)
            gcf_title = 'GCF(id={}), scores <strong>{}</strong>, shared strains=<strong>{}</strong>'.format(gcf.gcf_id, score, len(link.shared_strains))
            gcf_body = self.generate_gcf_info(gcf, link=link)
            span_badges = '&nbsp|&nbsp;'.join(['<span class="badge badge-info">{}</span>'.format(m.name) for m in link.methods])
            body += TMPL.format(hdr_id=gcf_hdr_id, hdr_color='ffe0b5', btn_target=gcf_body_id, btn_text=gcf_title, span_badges=span_badges,
                                body_id=gcf_body_id, body_parent='accordion_spec_{}'.format(pgindex), body_body=gcf_body)

            j += 1

        body += '</div>'
        return body

    def update_results_gcf_spec(self):
        """
        Generates scoring output for mode BGC/GCF -> Spectra. 

        This will be of the form:
            GCF 1
                <collapsible list of spectra linked to GCF 1>
            GCF 2
                <collapsible list of spectra linked to GCF 2>
            ...
        """
        pgindex = 0
        unsorted = []

        # iterate over every GCF and its associated set of linked objects + scores
        for gcf, link_data in self.score_helper.scoring_results.links.items():
            hdr_id = 'gcf_result_header_{}'.format(pgindex)
            body_id = 'gcf_result_body_{}'.format(pgindex)

            title = '{} spectra linked to GCF(id={})'.format(len(link_data), gcf.gcf_id)

            # first part of the body content is basic info about the GCF itself
            body = self.generate_gcf_info(gcf)
            body += '<hr/><h5>Linked objects:</h5>'
            # the second part is the list of linked spectra for the GCF, which is
            # generated by calling this method
            body += self.generate_gcf_spec_result(pgindex, gcf, link_data)

            # finally, store the complete generated HTML string in a list along with
            # the number of links this GCF has, so we can sort the eventual list
            # by number of links 
            unsorted.append((len(link_data), TMPL.format(hdr_id=hdr_id, hdr_color='adeaad', 
                                                           btn_target=body_id, btn_text=title, span_badges='',
                                                           body_id=body_id, body_parent='accordionResults', 
                                                           body_body=body)))
            pgindex += 1

        content = '<h4>GCFs={}, linked spectra={}</h4>'.format(len(self.score_helper.scoring_results), len(link_data))

        for _, text in sorted(unsorted, key=lambda x: -x[0]):
            content += text

        return content

    def update_results_gcf_molfam(self):
        return '' # TODO

    def update_results_spec_gcf(self):
        pgindex = 0
        unsorted = []
        all_gcfs = set()

        # for each original input object
        for spec, link_data in self.score_helper.scoring_results.links.items():
            hdr_id = 'spec_result_header_{}'.format(pgindex)
            body_id = 'spec_result_body_{}'.format(pgindex)

            title = '{} GCFs linked to {}'.format(len(link_data), spec)

            all_gcfs.update(link_data.keys())

            body = ''
            body += self.generate_spec_info(spec)
            
            # set up the chemdoodle plot so it appears when the entry is expanded 
            btn_id = 'spec_btn_{}'.format(pgindex)
            plot_id = 'spec_plot_{}'.format(pgindex)
            body += '<center><canvas id="{}"></canvas></center>'.format(plot_id)
            # note annoying escaping required here, TODO better way of doing this?
            onclick = 'setupPlot(\'{}\', \'{}\', \'{}\');'.format(btn_id, plot_id, spec.to_jcamp_str())
    
            body += '<hr/><h5>Linked objects:</h5>'
            body += self.generate_spec_gcf_result(pgindex, spec, link_data)

            unsorted.append((len(link_data), TMPL_ON_CLICK.format(hdr_id=hdr_id, hdr_color='b4c4e8', 
                                                         btn_target=body_id, btn_onclick=onclick, btn_id=btn_id, 
                                                         btn_text=title, span_badges='', body_id=body_id, 
                                                         body_parent='accordionResults', body_body=body)))

            pgindex += 1
        
        content = '<h4>Spectra={}, linked GCFs={}</h4>'.format(len(self.score_helper.scoring_results), len(all_gcfs))

        for _, text in sorted(unsorted, key=lambda x: -x[0]):
            content += text

        return content

    def update_results_molfam_gcf(self):
        return '' # TODO

    def update_results(self, div):
        """
        Called after scoring is completed. Generates HTML to display results and 
        assigns it to the content of a bokeh div to insert into the DOM
        """
        content = ''
        if self.score_helper.mode == SCO_MODE_BGC_SPEC or self.score_helper.mode == SCO_MODE_GCF_SPEC:
            # both of these have the same format
            self.debug_log('update_results_gcf_spec')
            content = self.update_results_gcf_spec()
        elif self.score_helper.mode == SCO_MODE_BGC_MOLFAM or self.score_helper.mode == SCO_MODE_GCF_MOLFAM:
            # as do these
            self.debug_log('update_results_gcf_molfam')
            content = self.update_results_gcf_molfam()
        elif self.score_helper.mode == SCO_MODE_SPEC_GCF:
            self.debug_log('update_results_spec_gcf')
            content = self.update_results_spec_gcf()
        elif self.score_helper.mode == SCO_MODE_MOLFAM_GCF:
            self.debug_log('update_results_molfam_gcf')
            content = self.update_results_molfam_gcf()
        else:
            self.debug_log('Unknown scoring mode selected! {}'.format(self.score_helper.mode))

        div.text = content

    def display_links(self, scoring_results, mode, div, include_only=None):
        # 1. need to filter out linked objects that don't appear in the tables 
        # in their current state
        if include_only is not None:
            scoring_results.filter_targets(lambda x: x in include_only)

        # TODO should probably remove this, or at least make it optional
        # 2. if the linked objects are GCFs, need to filter out those results for
        # which the GCFs are not currently "available" due to possible use of a
        # limited projection (this is tied to old TSNE plot projection stuff but
        # keeping it here in case still useful)
        # if mode == SCO_MODE_SPEC_GCF or mode == SCO_MODE_MOLFAM_GCF:
        #     scoring_results.filter_targets(lambda x: x in self.nh.available_gcfs[self.bgc_tsne_id])

        # 3. filter out any results with no shared strains, if that option is enabled
        if self.filter_no_shared_strains:
            scoring_results.filter_no_shared_strains()

        if len(scoring_results) > 0:
            self.debug_log('set_results')
            self.score_helper.set_results(scoring_results)
            self.update_alert('{} objects with links found'.format(len(self.score_helper.scoring_results)), 'primary')
        else:
            self.score_helper.clear()
            self.debug_log('No links found!')
            self.update_alert('No links found for last selection', 'danger')

        self.debug_log('update_results')
        self.update_results(div)

    def metcalf_standardised_callback(self, attr, old, new):
        obj = self.scoring_objects['metcalf']
        # don't attempt to set the same value it already has 
        # (this is using 0 = False, 1 = True)
        if obj.standardised == new:
            return
        obj.standardised = (new == 1)

    def metcalf_cutoff_callback(self, attr, old, new):
        self.scoring_objects['metcalf'].cutoff = new

    def rosetta_bgc_cutoff_callback(self, attr, old, new):
        self.scoring_objects['rosetta'].bgc_score_cutoff = new

    def rosetta_spec_cutoff_callback(self, attr, old, new):
        self.scoring_objects['rosetta'].spec_score_cutoff = new
        
    def scoring_mode_toggle_callback(self, attr, old, new):
        # 0 = AND, 1 = OR
        self.score_helper.method_and_mode = (new == 0)
        print('AND mode? {}'.format(self.score_helper.method_and_mode))

    def scoring_method_callback(self, attr, old, new):
        # update visibility of buttons
        for i in range(len(self.scoring_methods)):
            self.scoring_method_buttons[self.scoring_methods[i]].visible = i in new

    def get_scoring_mode_text(self):
        return SCO_MODE_NAMES[self.score_helper.mode]

    def sco_mode_changed(self):
        self.update_alert('Scoring mode is now <strong>{}</strong>'.format(self.get_scoring_mode_text()))

    def update_alert(self, msg, alert_class='primary'):
        self.alert_div.text = '<div class="alert alert-{}" role="alert">{}</div>'.format(alert_class, msg)

    def clear_alert(self):
        self.alert_div.text = ''

    def debug_log(self, msg, append=True):
        if append:
            self.debug_div.text += '{}<br/>'.format(msg)
        else:
            self.debug_div.text = '{}</br>'.format(msg)
        print('debug_log: {}'.format(msg))

    def search_type_callback(self, attr, old, new):
        print('Search type is now: {}={}'.format(new, SEARCH_OPTIONS.index(new)))

    def generate_bgc_info(self, bgc):
        bgc_body = '<ul>'
        for attr in ['strain', 'name', 'bigscape_class', 'product_prediction', 'description']:
            bgc_body += '<li><strong>{}</strong>: {}</li>'.format(attr, getattr(bgc, attr))
        bgc_body += '</ul>'
        return bgc_body

    def generate_gcf_info(self, gcf, link=None):
        shared_strains = []
        rosetta_obj = None
        if link is not None:
            shared_strains = link.shared_strains
            tmp = list(filter(lambda m: m.name == 'rosetta', link.methods))
            if len(tmp) > 0:
                rosetta_obj = tmp[0]

        gcf_body = '<ul class="nav nav-tabs" id="gcf_{}_tabs" role="tablist">'.format(gcf.id)
        gcf_body += '<li class="nav-item"><a class="nav-link active" id="gcf_{}_main_tab" data-toggle="tab" href="#gcf_{}_main" role="tab">Main</a></li>'.format(gcf.id, gcf.id)
        gcf_body += '<li class="nav-item"><a class="nav-link" id="gcf_{}_bgcs_tab" data-toggle="tab" href="#gcf_{}_bgcs" role="tab">BGCs</a></li>'.format(gcf.id, gcf.id)
        if rosetta_obj is not None:
            gcf_body += '<li class="nav-item"><a class="nav-link" id="gcf_{}_rosetta_tab" data-toggle="tab" href="#gcf_{}_rosetta" role="tab">Rosetta hits</a></li>'.format(gcf.id, gcf.id)
        gcf_body += '</ul><div class="tab-content" id="gcf_{}_tab_content">'.format(gcf.id)

        # start of main tab content
        gcf_body += '<div class="tab-pane show active" id="gcf_{}_main" role="tabpanel">'.format(gcf.id)
        gcf_body += '<ul>'
        for attr in ['id', 'gcf_id', 'bigscape_class']:
            gcf_body += '<li><strong>{}</strong>: {}</li>'.format(attr, getattr(gcf, attr))

        # add strain information
        if shared_strains is not None and len(shared_strains) > 0:
            gcf_body += '<li><strong>strains (total={}, shared={})</strong>: '.format(len(gcf.bgcs), len(shared_strains))

            non_shared = [s.strain for s in gcf.bgcs if s.strain not in shared_strains]

            for s in shared_strains:
                gcf_body += '<span style="background-color: #AAFFAA">{}</span>, '.format(s.id)
            for s in non_shared:
                gcf_body += '<span style="background-color: #DDDDDD">{}</span>, '.format(s.id)

            gcf_body += '</li>'
        else:
            gcf_body += '<li><strong>strains (total={}, shared=0)</strong>: '.format(len(gcf.bgcs))

            for s in gcf.bgcs:
                gcf_body += '<span>{}</span>, '.format(s.strain.id)

            gcf_body += '</li>'

        gcf_body += '</ul>'
        gcf_body += '</div>'

        # BGC tab content
        fields = ['name', 'strain', 'description', 'bigscape_classes', 'product_prediction']
        gcf_body += '<div class="tab-pane" id="gcf_{}_bgcs" role="tabpanel">'.format(gcf.id)
        gcf_body += 'Type to filter: <input type="text" onkeyup="onChangeBGCTable(\'#gcf_{}_bgc_table\', \'#gcf_{}_bgc_search\')" id="gcf_{}_bgc_search">'.format(gcf.id, gcf.id, gcf.id)
        gcf_body += '<table class="table table-responsive table-striped" id="gcf_{}_bgc_table"><thead><tr>'.format(gcf.id)
        gcf_body += ''.join('<th scope="col">{}</th>'.format(x) for x in fields)
        gcf_body += '</thead><tbody>'
        for bgc in gcf.bgcs:
            gcf_body += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(*list(getattr(bgc, x) for x in fields))
        gcf_body += '</tbody></table>'
        gcf_body += '</div>'

        # rosetta tab content, if available
        if rosetta_obj is not None:
            rosetta_fields = ['BGC name', 'BGC strain', 'MiBIG ID', 'BGC score', 'Spectrum ID', 'GNPS ID', 'Spectrum score']
            gcf_body += '<div class="tab-pane" id="gcf_{}_rosetta" role="tabpanel">'.format(gcf.id)
            gcf_body += 'Type to filter: <input type="text" onkeyup="onChangeRosettaTable(\'#gcf_{}_rosetta_table\', \'#gcf_{}_rosetta_search\')" id="gcf_{}_rosetta_search">'.format(gcf.id, gcf.id, gcf.id)
            gcf_body += '<table class="table table-responsive table-striped" id="gcf_{}_rosetta_table"><thead><tr>'.format(gcf.id)
            gcf_body += ''.join('<th scope="col">{}</th>'.format(x) for x in rosetta_fields)
            gcf_body += '</thead><tbody>'
            for hit in link.data(rosetta_obj):
                bgc_url = '<a href="https://mibig.secondarymetabolites.org/repository/{}" target="_blank">{}</a>'.format(hit.mibig_id, hit.mibig_id)
                gcf_body += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.3f}</td>'.format(hit.bgc.name, hit.bgc.strain, bgc_url, hit.bgc_match_score)
                spec_url = '<a href="{}" target="_blank">{}</a>'.format(gnps_url(hit.gnps_id, 'spectrum'), hit.gnps_id)
                gcf_body += '<td>{}</td><td>{}</td><td>{:.3f}</td></tr>'.format(hit.spec.spectrum_id, spec_url, hit.spec_match_score)
            gcf_body += '</tbody></table>'
            gcf_body += '</div>'

        gcf_body += '</div>'

        return gcf_body

    def generate_spec_info(self, spec, link=None):
        rosetta_obj = None
        if link is not None:
            tmp = list(filter(lambda m: m.name == 'rosetta', link.methods))
            if len(tmp) > 0:
                rosetta_obj = tmp[0]

        spec_body = '<ul class="nav nav-tabs" id="spec_{}_tabs" role="tablist">'.format(spec.id)
        spec_body += '<li class="nav-item"><a class="nav-link active" id="spec_{}_main_tab" data-toggle="tab" href="#spec_{}_main" role="tab">Main</a></li>'.format(spec.id, spec.id)
        if len(spec.annotations) > 0:
            spec_body += '<li class="nav-item"><a class="nav-link" id="spec_{}_anno_tab" data-toggle="tab" href="#spec_{}_anno" role="tab">Annotations</a></li>'.format(spec.id, spec.id)
        if rosetta_obj is not None:
            spec_body += '<li class="nav-item"><a class="nav-link" id="spec_{}_rosetta_tab" data-toggle="tab" href="#spec_{}_rosetta" role="tab">Rosetta hits</a></li>'.format(spec.id, spec.id)
        spec_body += '</ul><div class="tab-content" id="spec_{}_tab_content">'.format(spec.id)

        # main tab content
        spec_body += '<div class="tab-pane show active" id="spec_{}_main" role="tabpanel">'.format(spec.id)
        spec_body += '<ul>'
        for attr in ['id', 'spectrum_id', 'family', 'rt', 'total_ms2_intensity', 'max_ms2_intensity', 'n_peaks', 'precursor_mz', 'parent_mz']:
            spec_body += '<li><strong>{}</strong>: {}</li>'.format(attr, getattr(spec, attr))

        # TODO possibly going to need updated later
        if 'ATTRIBUTE_SampleType' in spec.metadata:
            spec_body += '<li><strong>SampleType</strong>: {}</li>'.format(', '.join(spec.metadata['ATTRIBUTE_SampleType']))

        # add strain information (possibly empty)
        if link is not None and len(link.shared_strains) > 0:
            spec_body += '<li><strong>strains (total={}, shared={})</strong>: '.format(len(spec.strains), len(link.shared_strains))

            non_shared = [s for s in spec.strains if s not in link.shared_strains]

            for s in link.shared_strains:
                spec_body += '<span style="background-color: #AAFFAA">{} ({})</span>, '.format(s.id, spec.get_growth_medium(s))
            for s in non_shared:
                spec_body += '<span style="background-color: #DDDDDD">{} ({})</span>, '.format(s.id, spec.get_growth_medium(s))

            spec_body += '</li>'
        else:
            spec_body += '<li><strong>strains (total={}, shared=0)</strong>: '.format(len(spec.strains))

            for s in spec.strains:
                spec_body += '<span>{} ({})</span>, '.format(s.id, spec.get_growth_medium(s))

            spec_body += '</li>'
        spec_body += '</ul></div>'

        # annotations tab
        if len(spec.annotations) > 0:
            spec_body += '<div class="tab-pane" id="spec_{}_anno" role="tabpanel">'.format(spec.id)

            # keys of spec.annotations indicate the source file: "gnps" or an actual filename
            for anno_src, anno_list in spec.annotations.items():
                spec_body += '<strong>{} annotations:</strong><ul>'.format(anno_src)
                for item in anno_list:
                    spec_body += '<ul>'
                    content = []
                    for k, v in item.items():
                        if k.endswith('_url'):
                            v = '<a href="{}" target="_blank">{}</a>'.format(v, k)
                        elif v.startswith('CCMSLIB'):
                            v = '<a href="{}" target="_blank">{}</a>'.format(gnps_url(v), v)
                        content.append('<strong><span class="annotation">{}</span></strong> = {}'.format(k, v))
                    spec_body += '<li>' + ', '.join(content) + '</li>'
                    spec_body += '</ul>'
                spec_body += '</ul>'
                spec_body += '<hr width="80%"/>'

            spec_body += '</div>'

        # rosetta tab content, if available
        if rosetta_obj is not None:
            rosetta_fields = ['BGC name', 'BGC strain', 'MiBIG ID', 'BGC score', 'Spectrum ID', 'GNPS ID', 'Spectrum score']
            spec_body += '<div class="tab-pane" id="spec_{}_rosetta" role="tabpanel">'.format(spec.id)
            spec_body += 'Type to filter: <input type="text" onkeyup="onChangeRosettaTable(\'#spec_{}_rosetta_table\', \'#spec_{}_rosetta_search\')" id="spec_{}_rosetta_search">'.format(spec.id, spec.id, spec.id)
            spec_body += '<table class="table table-responsive table-striped" id="spec_{}_rosetta_table"><thead><tr>'.format(spec.id)
            spec_body += ''.join('<th scope="col">{}</th>'.format(x) for x in rosetta_fields)
            spec_body += '</thead><tbody>'
            for hit in link.data(rosetta_obj):
                bgc_url = '<a href="https://mibig.secondarymetabolites.org/repository/{}" target="_blank">{}</a>'.format(hit.mibig_id, hit.mibig_id)
                spec_body += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.3f}</td>'.format(hit.bgc.name, hit.bgc.strain, bgc_url, hit.bgc_match_score)
                spec_url = '<a href="{}" target="_blank">{}</a>'.format(gnps_url(hit.gnps_id, 'spectrum'), hit.gnps_id)
                spec_body += '<td>{}</td><td>{}</td><td>{:.3f}</td></tr>'.format(hit.spec.spectrum_id, spec_url, hit.spec_match_score)
            spec_body += '</tbody></table>'
            spec_body += '</div>'

        spec_body += '</div>'
        return spec_body

    def generate_molfam_info(self, molfam, shared_strains=[]):
        pass # TODO

    def generate_search_output_bgc(self, bgcs):
        body = ''
        for i, bgc in enumerate(bgcs):
            bgc_hdr_id = 'bgc_search_header_{}'.format(i)
            bgc_body_id = 'bgc_search_body_{}'.format(i)
            bgc_title = '{}(name={})'.format(bgc.__class__.__name__, bgc.name)
            bgc_body = self.generate_bgc_info(bgc)
            body += TMPL_SEARCH.format(hdr_id=bgc_hdr_id, hdr_color='dddddd', btn_target=bgc_body_id, btn_text=bgc_title, 
                                result_index=str(i), body_id=bgc_body_id, body_parent='accordionSearch', body_body=bgc_body)

        
        return body

    def generate_search_output_gcf(self, gcfs):
        body = ''
        for i, gcf in enumerate(gcfs):
            gcf_hdr_id = 'gcf_search_header_{}'.format(i)
            gcf_body_id = 'gcf_search_body_{}'.format(i)
            gcf_title = 'GCF(gcf_id={}, strains={})'.format(gcf.gcf_id, len(gcf.bgcs))
            gcf_body = self.generate_gcf_info(gcf)
            body += TMPL_SEARCH.format(hdr_id=gcf_hdr_id, hdr_color='dddddd', btn_target=gcf_body_id, btn_text=gcf_title, 
                                result_index=str(i), body_id=gcf_body_id, body_parent='accordionSearch', body_body=gcf_body)

        return body


    def generate_search_output_spec(self, spectra):
        body = ''
        for i, spec in enumerate(spectra):
            spec_hdr_id = 'spec_search_header_{}'.format(i)
            spec_body_id = 'spec_search_body_{}'.format(i)
            spec_title = 'Spectrum(parent_mass={:.4f}, id={}, strains={})'.format(spec.parent_mz, spec.spectrum_id, len(spec.strains))

            spec_body = self.generate_spec_info(spec)

            # set up the chemdoodle plot so it appears when the entry is expanded 
            spec_btn_id = 'spec_search_btn_{}'.format(i)
            spec_plot_id = 'spec_search_plot_{}'.format(i)
            spec_body += '<center><canvas id="{}"></canvas></center>'.format(spec_plot_id)

            # note annoying escaping required here, TODO better way of doing this?
            spec_onclick = 'setupPlot(\'{}\', \'{}\', \'{}\');'.format(spec_btn_id, spec_plot_id, spec.to_jcamp_str())

            body += TMPL_SEARCH_ON_CLICK.format(hdr_id=spec_hdr_id, hdr_color='dddddd', btn_target=spec_body_id, btn_onclick=spec_onclick, btn_id=spec_btn_id, 
                                         result_index=str(i), btn_text=spec_title, body_id=spec_body_id, body_parent='accordionSearch', body_body=spec_body)

        return body

    def generate_search_output_molfam(self, molfams):
        # TODO
        return '<strong>Not working yet :(</strong>'

    def generate_search_output(self, results):
        hdr = '<h4>{} results found</h4>'.format(len(results))
        if len(results) == 0:
            self.search_div_header.text = hdr
            self.search_div_results.text = ''
            return

        hdr = '<h4><button type="button" class="btn btn-info" onclick="toggleSearchSelect()">{} results found, click to select/deselect all</button></h4>'.format(len(results))
        self.search_div_header.text = hdr

        html = ''
        if isinstance(results[0], BGC):
            html += self.generate_search_output_bgc(results)
        elif isinstance(results[0], GCF):
            html += self.generate_search_output_gcf(results)
        elif isinstance(results[0], Spectrum):
            html += self.generate_search_output_spec(results)
        elif isinstance(results[0], MolecularFamily):
            html += self.generate_search_output_molfam(results)
            
        self.search_div_results.text = html

    def search_button_callback(self):
        self.generate_search_output(self.searcher.search(self.search_type.value, self.search_input.value, len(self.search_regex.active) == 1))

    def search_score_button_callback(self):
        # check if the user responded yes/no to query (the JS callback for this button 
        # is triggered before the python callback)
        if not self.confirm_response.data['resp'][0]:
            print('Bailing out due to JS confirmation')
            return

        # retrieve the set of selected indices from the Javascript callback and create subset of results
        selected_indices = list(map(int, self.confirm_response.data['selected'][0]))
        results = [self.searcher.results[i] for i in selected_indices]
        if len(results) == 0:
            self.update_alert('No search results!')
            return

        # may need to change scoring mode here
        # if results are BGCs, switch to genomics mode and select the BGC mode
        #
        # TODO: need to fix this for tables scoring modes
        if isinstance(results[0], BGC):
            print('scoring based on bgc results')
            print('current scoring mode gen? {}'.format(self.score_helper.from_genomics))
            # now select bgc mode
            self.score_helper.update_genomics(GENOMICS_MODE_BGC)

            indices = []
            for bgc in results:
                for i in range(len(self.table_session_data.bgc_ds.data['bgc_pk'])):
                    if bgc.id == self.table_session_data.bgc_ds.data['bgc_pk'][i]:
                        indices.append(i)
                        break

            self.table_session_data.bgc_ds.selected.indices = indices
            self.tables_score_gen_callback()

        # if results are GCFs, switch to genomics mode and select the GCF mode
        elif isinstance(results[0], GCF):
            print('scoring based on gcf results')
            print('current scoring mode gen? {}'.format(self.score_helper.from_genomics))
            # now select bgc mode
            self.score_helper.update_genomics(GENOMICS_MODE_BGC)

            # update selection to trigger scoring, by building a list
            # of all the BGC indices corresponding to these GCFs
            bgcs = set()
            for gcf in results:
                bgcs.update(gcf.bgcs)

            indices = []
            for bgc in bgcs:
                for i in range(len(self.table_session_data.bgc_ds.data['bgc_pk'])):
                    if bgc.id == self.table_session_data.bgc_ds.data['bgc_pk'][i]:
                        indices.append(i)
                        break

            self.table_session_data.bgc_ds.selected.indices = indices
            self.tables_score_gen_callback()
        elif isinstance(results[0], Spectrum):
            print('scoring based on spectra results')
            print('current scoring mode met? {}'.format(not self.score_helper.from_genomics))

            # now select spectrum mode
            self.score_helper.update_metabolomics(METABOLOMICS_MODE_SPEC)

            # need to find the spectra in the tables datasource and add to selected list
            indices = []
            for spec in results:
                for i in range(len(self.table_session_data.spec_ds.data['spectrum_id'])):
                    if spec.spectrum_id == self.table_session_data.spec_ds.data['spectrum_id'][i]:
                        indices.append(i)
                        break

            self.table_session_data.spec_ds.selected.indices = indices
            self.tables_score_met_callback()

    def gnps_params_select_callback(self, attr, old, new):
        self.gnps_params_value.text = '<strong>{}</strong> = {}'.format(new, self.nh.nplinker.gnps_params[new])

    def current_scoring_method_names(self):
        print('NAMES', [self.scoring_methods[x] for x in self.scoring_method_checkboxes.active])
        return [self.scoring_methods[x] for x in self.scoring_method_checkboxes.active]

    def current_scoring_methods(self):
        print('METHODS', [self.scoring_objects[x] for x in self.current_scoring_method_names()])
        return [self.scoring_objects[x] for x in self.current_scoring_method_names()]

    def tables_score_met_callback(self):
        """
        Handles clicks on the "Show scores for selected spectra" button (tables mode)
        """
        print('Currently selected spectra: {}'.format(len(self.table_session_data.spec_ds.data['spectrum_id'])))
        print('Actual selections: {}'.format(len(self.table_session_data.spec_ds.selected.indices)))

        if len(self.table_session_data.spec_ds.selected.indices) == 0 and len(self.table_session_data.molfam_ds.selected.indices) == 0:
            self.hidden_alert_thing.text = 'Error: no spectra have been selected!\n\n' +\
                                            'Running scoring on {} objects may take a long time!\n'.format(len(self.table_session_data.spec_ds.data['spectrum_id'])) +\
                                             'Please make some selections and then try again.'
            # this resets the text without triggering the alert() again, so that NEXT time this happens
            # the text-changed event will be fired by the above line. yes it's a mess
            self.hidden_alert_thing.text = ''
            return

        # need to build a list of the NPLinker Spectrum objects corresponding to
        # the set of objects currently listed in the table
        # but first need to decide which set of spectra to use:
        # - the user might have selected a MolFam, in which case the spec table will have
        #   filtered down to show only the spectra in that family, BUT none of them will
        #   actually be selected, so just take all visible entries
        # - alternatively they might have selected a spectrum directly, in which case they
        #   will all still be visible but we only want the selected indices
        selected_spectra = []
        if len(self.table_session_data.spec_ds.selected.indices) > 0:
            # assume there's been an explicit selection and use those objects
            for index in self.table_session_data.spec_ds.selected.indices:
                sid = self.table_session_data.spec_ds.data['spec_pk'][index] # the internal nplinker ID
                if sid == NA_ID:
                    continue
                selected_spectra.append(self.nh.nplinker.spectra[int(sid)])
        else:
            for sid in self.table_session_data.spec_ds.data['spec_pk']: # internal nplinker IDs
                if sid == NA_ID:
                    continue
                spec = self.nh.nplinker.spectra[int(sid)] 
                selected_spectra.append(spec)

        # retrieve the currently active scoring method(s)
        current_methods = self.current_scoring_methods()
        
        # call get_links, which returns the subset of input objects which have links
        results = self.nh.nplinker.get_links(selected_spectra, current_methods, and_mode=self.score_helper.method_and_mode)

        # TODO quick attempt to display scores in the tables (for metcalf only for now)
        # TODO does it make sense to do this with other methods? might not have a single score...
        metcalf_obj = self.scoring_objects['metcalf']
        table_metcalf_scores = ['-' for gcf_id in self.table_session_data.gcf_ds.data['gcf_id']]
        # check if metcalf scoring is enabled at the moment
        if metcalf_obj in current_methods:
            # create a dict mapping GCFs in the results to their scores
            # NOTE: rosetta scoring currently outputs BGCs not GCFs, which breaks
            # this assumption that GCFs are present. currently working around by 
            # setting bgc_to_gcf=True on the scoring object...
            gcf_scores = {}
            for spec, result in results.links.items():
                for gcf, link_data in result.items():
                    if metcalf_obj in link_data.methods:
                        metcalf_score = link_data.data(metcalf_obj)
                        gcf_scores[gcf] = metcalf_score

            # get a list of the GCFs currently visible in the table (note that
            # these are the actual GCF IDs, not the internal nplinker IDs)
            visible_gcfs = list(map(int, self.table_session_data.gcf_ds.data['gcf_id']))

            # update the data source with the scores
            table_metcalf_scores = []
            for gcf_id in visible_gcfs:
                gcf_obj = self.nh.nplinker.lookup_gcf(gcf_id)
                if gcf_obj not in gcf_scores:
                    table_metcalf_scores.append('-')
                else:
                    table_metcalf_scores.append(gcf_scores[gcf_obj])

        self.table_session_data.gcf_ds.data['metcalf_score'] = table_metcalf_scores

        # the next step, before actually displaying the links, is to make sure
        # the webapp scoring mode is set to Spec->GCF
        self.score_helper.update_metabolomics(METABOLOMICS_MODE_SPEC)
        self.sco_mode_changed()

        # the ScoringHelper object expects to have this list of objects available so need to set that
        self.score_helper.spectra = selected_spectra

        # get the list of visible GCFs so we can filter out any others from the results
        include_only = set([self.nh.nplinker.gcfs[int(gcf_id)] for gcf_id in self.table_session_data.gcf_ds.data['gcf_pk'] if gcf_id != '-'])

        # FINALLY, display the link information
        self.display_links(results, mode=SCO_MODE_SPEC_GCF, div=self.results_div, include_only=include_only)
        
    def tables_score_gen_callback(self):
        """
        Handles clicks on the "Show scores for selected GCFs" button (tables mode)
        """

        print('Currently selected bgcs: {}'.format(len(self.table_session_data.bgc_ds.data['bgc_pk'])))
        print('Actual selections: {}'.format(len(self.table_session_data.bgc_ds.selected.indices)))
        ignore_hybrid_bgcs = len(self.hybrid_scoring_control_checkbox.active) == 1 

        if len(self.table_session_data.bgc_ds.selected.indices) == 0 and len(self.table_session_data.gcf_ds.selected.indices) == 0:
            self.hidden_alert_thing.text = 'Error: no BGCs have been selected!\n\n' +\
                                            'Running scoring on {} objects may take a long time!\n'.format(len(self.table_session_data.bgc_ds.data['bgc_pk'])) +\
                                             'Please make some selections and then try again.'
            # this resets the text without triggering the alert() again, so that NEXT time this happens
            # the text-changed event will be fired by the above line. yes it's a mess
            self.hidden_alert_thing.text = ''
            return

        # need to build a list of the NPLinker GCF objects corresponding to
        # the set of BGC objects currently listed in the table
        # but first need to decide which set of objects to use:
        # - the user might have selected a GCF, in which case the BGC table will have
        #   filtered down to show only the BGCs in that family, BUT none of them will
        #   actually be selected, so just take all visible entries
        # - alternatively they might have selected a BGC directly, in which case they
        #   will all still be visible but we only want the selected indices
        selected_gcfs = []
        selected_bgcs = []
        gcfs = set() # to remove dups 

        if len(self.table_session_data.bgc_ds.selected.indices) > 0:
            # assume there's been an explicit selection and use those objects
            # (i.e. the user has selected BGCs directly instead of filtering by 
            # picking a GCF first)
            for index in self.table_session_data.bgc_ds.selected.indices:
                bgcid = self.table_session_data.bgc_ds.data['bgc_pk'][index] # internal nplinker ID
                if bgcid == NA_ID:
                    continue
                bgc = self.nh.nplinker.bgcs[int(bgcid)]
                selected_bgcs.append(bgc)
                # TODO should probably add a UI widget to allow for control over
                # hybrid BGCs. e.g. if you select a BGC that has 3 GCF parents
                # with different BiG-SCAPE classes, it should be possible to say
                # "only include the NRPS one"... 
                gcfs.update(bgc.parents)
        else:
            print(self.table_session_data.gcf_ds.selected.indices)
            # get the selected GCFs
            tmp_gcfs = []
            for index in self.table_session_data.gcf_ds.selected.indices:
                gcfid = self.table_session_data.gcf_ds.data['gcf_pk'][index] # internal nplinker ID
                tmp_gcfs.append(self.nh.nplinker.gcfs[gcfid])

            for bgcid in self.table_session_data.bgc_ds.data['bgc_pk']: # internal nplinker IDs
                if bgcid == NA_ID:
                    continue
                bgc = self.nh.nplinker.bgcs[int(bgcid)]
                selected_bgcs.append(bgc)
                if ignore_hybrid_bgcs and bgc.is_hybrid:
                    # if ignoring hybrid BGCs, only want to include parent GCFs with 
                    # matching BiG-SCAPE classes, so only add those that appear in 
                    # the set above (that we know were explicitly selected by the user)
                    for parent in bgc.parents:
                        if parent in tmp_gcfs:
                            gcfs.add(parent)
                else:
                    # just add all parents
                    gcfs.update(bgc.parents)

        selected_gcfs = list(gcfs)

        # retrieve the currently active scoring method (for now this will always be metcalf)
        current_methods = self.current_scoring_methods()
        
        # call get_links, which returns the subset of input objects which have links
        results = self.nh.nplinker.get_links(selected_gcfs, current_methods, and_mode=self.score_helper.method_and_mode)

        # TODO quick attempt to display scores in the tables (for metcalf only for now)
        # TODO does it make sense to do this with other methods? might not have a single score...
        metcalf_obj = self.scoring_objects['metcalf']
        table_metcalf_scores = ['-' for spec_id in self.table_session_data.spec_ds.data['metcalf_score']]
        # check if metcalf scoring is enabled at the moment
        if metcalf_obj in current_methods:
            # create a dict mapping spectra in the results to their scores
            # NOTE: rosetta scoring currently outputs BGCs not GCFs, which breaks
            # this assumption that GCFs are present. currently working around by 
            # setting bgc_to_gcf=True on the scoring object...
            spec_scores = {}
            for gcf, result in results.links.items():
                for spec, link_data in result.items():
                    if metcalf_obj in link_data.methods:
                        metcalf_score = link_data.data(metcalf_obj)
                        spec_scores[spec] = metcalf_score

            # get a list of the spectra currently visible in the table (note that
            # these are the actual Spectrum IDs, not the internal nplinker IDs)
            visible_spectra = list(map(int, self.table_session_data.spec_ds.data['spectrum_id']))

            # update the data source with the scores
            table_metcalf_scores = []
            for spec_id in visible_spectra:
                spec_obj = self.nh.nplinker.lookup_spectrum(spec_id)
                if spec_obj not in spec_scores:
                    table_metcalf_scores.append('-')
                else:
                    table_metcalf_scores.append(spec_scores[spec_obj])

        self.table_session_data.spec_ds.data['metcalf_score'] = table_metcalf_scores

        # the next step, before actually displaying the links, is to make sure
        # the webapp scoring mode is set to GCF->Spec
        self.score_helper.update_genomics(GENOMICS_MODE_BGC)
        self.sco_mode_changed()

        # the ScoringHelper object expects to have this list of objects available so need to set that
        self.score_helper.bgcs = selected_bgcs
        self.score_helper.gcfs = selected_gcfs

        # get the list of visible spectra so we can filter out any others from the results
        # TODO MolFam mode...
        include_only = set([self.nh.nplinker.spectra[int(spec_id)] for spec_id in self.table_session_data.spec_ds.data['spec_pk'] if spec_id != '-'])

        # FINALLY, display the link information
        self.display_links(results, mode=SCO_MODE_BGC_SPEC, div=self.results_div, include_only=include_only)

    def tables_reset_callback(self):
        for dt in self.table_session_data.data_tables.values():
            dt.selectable = True

    def bokeh_layout(self):
        widgets = []
        self.results_div = Div(text=RESULTS_PLACEHOLDER, name='results_div', width_policy='max', height_policy='fit')
        widgets.append(self.results_div)

        # scoring objects and related stuff
        # TODO this should be drawn from nplinker config, not hardcoded here? or config file
        self.scoring_methods = ['metcalf', 'rosetta']
        self.scoring_objects = {m: self.nh.nplinker.scoring_method(m) for m in self.scoring_methods}

        # checkbox list of scoring methods to enable/disable them
        self.scoring_method_checkboxes = CheckboxGroup(labels=[m for m in self.scoring_methods], active=[0], name='scoring_method_checkboxes')
        self.scoring_method_checkboxes.on_change('active', self.scoring_method_callback)
        widgets.append(self.scoring_method_checkboxes)

        # buttons to show the modal dialogs for configuring each scoring method
        self.scoring_method_buttons = {}
        self.scoring_method_buttons['metcalf'] = Button(label='Configure Metcalf scoring', name='metcalf_scoring_button', button_type='primary')
        self.scoring_method_buttons['rosetta'] = Button(label='Configure rosetta scoring', name='rosetta_scoring_button', button_type='primary', visible=False)
        widgets.extend(x for x in self.scoring_method_buttons.values())

        # basic JS callbacks to show the modals that allow you to configure the enabled scoring methods
        show_modal_callback = """ $('#' + modal_name).modal(); """
        self.scoring_method_buttons['metcalf'].js_on_click(CustomJS(args=dict(modal_name='metcalfModal'), code=show_modal_callback))
        self.scoring_method_buttons['rosetta'].js_on_click(CustomJS(args=dict(modal_name='rosettaModal'), code=show_modal_callback))

        # metcalf config objects
        self.metcalf_standardised = RadioGroup(labels=['Basic Metcalf scoring', 'Standardised Metcalf scoring'], name='metcalf_standardised', active=1 if self.scoring_objects['metcalf'].standardised else 0, sizing_mode='scale_width')
        self.metcalf_standardised.on_change('active', self.metcalf_standardised_callback)
        widgets.append(self.metcalf_standardised)
        self.metcalf_cutoff = Slider(start=0, end=10, value=int(self.scoring_objects['metcalf'].cutoff), step=1, title='cutoff value', name='metcalf_cutoff')
        self.metcalf_cutoff.on_change('value', self.metcalf_cutoff_callback)
        widgets.append(self.metcalf_cutoff)

        # same for rosetta
        self.rosetta_spec_cutoff = Slider(start=0, end=1, step=0.1, value=self.scoring_objects['rosetta'].bgc_score_cutoff, name='rosetta_spec_cutoff', title='spectral score cutoff value')
        self.rosetta_spec_cutoff.on_change('value', self.rosetta_spec_cutoff_callback)
        self.rosetta_bgc_cutoff = Slider(start=0, end=1, step=0.1, value=self.scoring_objects['rosetta'].bgc_score_cutoff, name='rosetta_bgc_cutoff', title='BGC score cutoff value')
        self.rosetta_bgc_cutoff.on_change('value', self.rosetta_bgc_cutoff_callback)
        widgets.append(self.rosetta_spec_cutoff)
        widgets.append(self.rosetta_bgc_cutoff)

        self.rosetta_dl_button = Button(label='Download all Rosetta hits for this dataset as a CSV file', name='rosetta_dl_button')
        rosetta_download_code = """
            const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })
            const filename = 'rosetta_hits.csv';

            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename)
            } else {
                const link = document.createElement('a')
                link.href = URL.createObjectURL(blob)
                link.download = filename
                link.target = '_blank'
                link.style.visibility = 'hidden'
                link.dispatchEvent(new MouseEvent('click'))
            }
        """

        # load rosetta hits (if available)
        rosetta_path = os.path.join(self.nh.nplinker.root_dir, 'rosetta', 'rosetta_hits.csv')
        if os.path.exists(rosetta_path):
            with open(rosetta_path, 'r') as f:
                rosetta_csv_text = f.read()
        self.rosetta_dl_button.js_on_click(CustomJS(args=dict(filetext=rosetta_csv_text), code=rosetta_download_code))
        widgets.append(self.rosetta_dl_button)

        # and/or mode selection
        self.scoring_mode_toggle = RadioButtonGroup(labels=['Show results from all methods (AND mode)', 'Show results from any method (OR mode)'], active=0, name='scoring_mode_toggle', sizing_mode='scale_width')
        self.scoring_mode_toggle.on_change('active', self.scoring_mode_toggle_callback)
        widgets.append(self.scoring_mode_toggle)

        metcalf_info_div = Div(name='metcalf_info', sizing_mode='scale_width', text='Tables contain objects with Metcalf score > {:.1f}'.format(self.nh.nplinker._loader.webapp_scoring_cutoff()))
        self.metcalf_info = wrap_popover(metcalf_info_div, create_popover(*TOOLTIP_TABLES_FILTERING), 'metcalf_info')
        widgets.append(self.metcalf_info)

        # for debug output etc 
        self.debug_div = Div(text="", name='debug_div', sizing_mode='scale_width')
        widgets.append(self.debug_div)

        # status updates/errors
        self.alert_div = Div(text="", name='alert_div')
        widgets.append(self.alert_div)

        # header text
        self.header_text = Div(text='<strong>NPLinker web app</strong> (loaded BGCs={}, GCFs={}, Spectra={}, MolFams={})'.format(len(self.nh.nplinker.bgcs), len(self.nh.nplinker.gcfs), len(self.nh.nplinker.spectra), len(self.nh.nplinker.molfams)), name='header_text', sizing_mode='scale_width')
        widgets.append(self.header_text)

        self.search_input = TextInput(value='', title='Search for:', name='search_input', sizing_mode='scale_width')
        # TODO this can be used to trigger a search too, but annoyingly it is triggered when
        # a) enter is pressed (which is fine) 
        #   OR
        # b) the widget is unfocused (which is a bit stupid because if the user is unfocusing to
        # go click the search button it will end up doing the search twice)
        # self.search_input.on_change('value', lambda a, o, n : self.search_button_callback())
        widgets.append(self.search_input)

        # TODO can use this to convert selections in the table to selections on the plot
        self.foo = ColumnDataSource(data={'data': [1, 2, 3, 4, 5]})
        # self.foo.on_change('data', lambda a, o, n: print('FOO', a, o, n))
        self.search_button = Button(name='search_button', label='Search', button_type='primary', sizing_mode='scale_width')
        # conveniently on_click seems to get called after js_on_click, so we can do stuff on the 
        # JS side in the CustomJS callback and then pick up the changes in the on_click handler
        # on the Python side

        self.search_button.on_click(self.search_button_callback)
        widgets.append(self.search_button)

        self.search_score_button = Button(name='search_score_button', label='Run scoring on selected results', button_type='warning', sizing_mode='scale_width')
        self.search_score_button.on_click(self.search_score_button_callback)
        self.confirm_response = ColumnDataSource(data={'resp': [None], 'selelected': [[]]})
        self.search_score_button.js_on_click(CustomJS(args=dict(resp=self.confirm_response), code="""
                // purpose of this callback is to ask the user to confirm scoring
                // is desired if the number of selected objects is > 25 (arbitrary threshold),
                // and also to extract the indices of the selected objects to return to the
                // python code (this code is executed before the python callback)

                // get the array of selected indices
                var sel_indices = $('input[class=search-input]:checked').map(function() { return this.id; }).get();
                if(sel_indices.length > 25) {
                    // use a JS confirm dialog to get a response
                    var answer = confirm(sel_indices.length + " objects selected, do you really want to run scoring? (may be slow!)");
                    // switch to scoring tab if confirmed to proceed
                    if(answer == true) {
                        $('#scoringTab').tab('show');
                    }
                    // return the answer in the data source
                    resp.data = {'resp': [answer], 'selected': [sel_indices] };
                } else if(sel_indices.length > 0) {
                    // in this case can just proceed directly to show scoring tab and return selection
                    resp.data = {'resp': [true], 'selected': [sel_indices] };
                    $('#scoringTab').tab('show');
                }
                resp.change.emit();
                """))
        widgets.append(self.search_score_button)

        self.search_div_results = Div(text='', sizing_mode='scale_both', name='search_div_results')
        widgets.append(self.search_div_results)

        self.search_div_header = Div(text='', sizing_mode='scale_height', name='search_div_header')
        widgets.append(self.search_div_header)

        self.search_type = Select(options=SEARCH_OPTIONS, value=SEARCH_OPTIONS[0], name='search_type', css_classes=['nolabel'], sizing_mode='scale_width')
        self.search_type.on_change('value', self.search_type_callback)
        widgets.append(self.search_type)

        self.search_regex = CheckboxGroup(labels=['Use regexes?'], active=[0], name='search_regex', sizing_mode='scale_width')
        widgets.append(self.search_regex)

        self.dataset_description_pre = PreText(text=self.nh.nplinker.dataset_description, sizing_mode='scale_both', name='dataset_description_pre')
        widgets.append(self.dataset_description_pre)

        gnps_params = self.nh.nplinker.gnps_params 
        defval = 'uuid' if 'uuid' in gnps_params else None
        self.gnps_params_select = Select(options=list(gnps_params.keys()), value=defval, name='gnps_params_select', css_classes=['nolabel'], sizing_mode='scale_height')
        self.gnps_params_select.on_change('value', self.gnps_params_select_callback)
        widgets.append(self.gnps_params_select)

        if defval is not None:
            deftext = '<strong>{}</strong> = {}'.format(defval, self.nh.nplinker.gnps_params[defval])
        else:
            deftext = ''
        self.gnps_params_value = Div(text=deftext, sizing_mode='scale_height', name='gnps_params_value')
        widgets.append(self.gnps_params_value)

        other_info = '<table class="table table-responsive table-striped" id="other_info_table"><thead><tr>'
        other_info += '<th scope="col">Name</th><th scope="col">Value</th></tr>'
        other_info += '</thead><tbody>'
        other_info += '<tr><td>{}</td><td>{}</td></tr>'.format('BiGSCAPE clustering cutoff', self.nh.nplinker.bigscape_cutoff)
        other_info += '</tbody></table>'
        self.other_info_div = Div(text=other_info, sizing_mode='scale_height', name='other_info_div')
        widgets.append(self.other_info_div)

        # tables widgets
        self.table_session_data.setup()
        widgets.extend(self.table_session_data.widgets)

        # callbacks for scoring buttons
        self.table_session_data.tables_score_met.on_click(self.tables_score_met_callback)
        self.table_session_data.tables_score_gen.on_click(self.tables_score_gen_callback)

        # and reset button
        self.table_session_data.tables_reset_button.on_click(self.tables_reset_callback)

        self.hidden_alert_thing = PreText(text='', name='hidden_alert', visible=False)
        self.hidden_alert_thing.js_on_change('text', CustomJS(code='if(cb_obj.text.length > 0) alert(cb_obj.text);', args={}))
        widgets.append(self.hidden_alert_thing)

        # hybrid scoring control
        self.hybrid_scoring_control_checkbox = CheckboxGroup(sizing_mode='scale_width', labels=['Treat hybrid BGCs as non-hybrid during scoring'], active=[0], name='hybrid_scoring_control_checkbox')
        self.hybrid_scoring_control_checkbox_wrapper = wrap_popover(self.hybrid_scoring_control_checkbox, create_popover(*TOOLTIP_HYBRID_SCORING), 'hybrid_scoring_control_checkbox', sizing_mode='scale_width')
        widgets.append(self.hybrid_scoring_control_checkbox_wrapper)

        print('layout done, adding widget roots')

        for w in widgets:
            curdoc().add_root(w)

        curdoc().title = 'nplinker webapp'
        print('bokeh_layout method complete')
        # curdoc().theme = 'dark_minimal'

# server_lifecycle.py adds a .nh attr to the current Document instance, use that
# to access the already-created NPLinker instance plus the TSNE data
nb = NPLinkerBokeh(curdoc().nh, curdoc().table_session_data)
nb.bokeh_layout()
nb.update_alert('Initialised OK!', 'success')
