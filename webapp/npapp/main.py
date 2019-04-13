from bokeh.models.widgets import RadioGroup, Slider, Div, CheckboxButtonGroup, Select, CheckboxGroup, Toggle, Button
from bokeh.layouts import row, column, widgetbox
from bokeh.models import CustomJS
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, HoverTool, TapTool, CDSView, IndexFilter
from bokeh.events import LODEnd, LODStart, Reset

from nplinker import Spectrum, GCF, MolecularFamily

# TOOLS="crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
TOOLS = "crosshair,pan,wheel_zoom,box_zoom,reset,save,box_select,tap"

PW, PH = 500, 500

POINT_FACTOR = 150.0
POINT_RATIO = 1 / POINT_FACTOR

ACTIVE_BG_FILL = '#ffffff'
INACTIVE_BG_FILL = '#eeeeee'
ACTIVE_BG_ALPHA = 1.0
INACTIVE_BG_ALPHA = 0.555

SCORING_GEN_TO_MET = 'TO <------ FROM'
SCORING_MET_TO_GEN = 'FROM ------> TO'

PLOT_TOGGLES = ['Enable alpha-blending', 'Enable colormaps', 'Show singleton families', 'Preserve colours when selecting', 'Only show results with shared strains']
PLOT_ALPHA, PLOT_CMAP, PLOT_SINGLETONS, PLOT_PRESERVE_COLOUR, PLOT_ONLY_SHARED_STRAINS = range(len(PLOT_TOGGLES))
PLOT_TOGGLES_ENUM = [PLOT_ALPHA, PLOT_CMAP, PLOT_SINGLETONS, PLOT_PRESERVE_COLOUR, PLOT_ONLY_SHARED_STRAINS]

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

class ScoringHelper(object):

    GENOMICS_LOOKUP  = {GENOMICS_MODE_BGC: [SCO_MODE_BGC_SPEC, SCO_MODE_BGC_MOLFAM],
                        GENOMICS_MODE_GCF: [SCO_MODE_GCF_SPEC, SCO_MODE_GCF_MOLFAM]}

    METABOLOMICS_LOOKUP = [SCO_MODE_SPEC_GCF, SCO_MODE_MOLFAM_GCF]

    def __init__(self, nh):
        self.genomics_mode = GENOMICS_MODE_BGC
        self.metabolomics_mode = METABOLOMICS_MODE_SPEC
        self.mode = SCO_MODE_BGC_SPEC
        self.from_genomics = True
        self.nh = nh
        self.clear()

    def clear(self):
        self.bgcs = []
        self.gcfs = []
        self.spectra = []
        self.molfams = []
        self.inputs = set()

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

    def _update(self):
        if self.from_genomics:
            self.mode = ScoringHelper.GENOMICS_LOOKUP[self.genomics_mode][self.metabolomics_mode]
        else:
            self.mode = ScoringHelper.METABOLOMICS_LOOKUP[self.metabolomics_mode]

    @property
    def mode_name(self):
        return SCO_MODE_NAMES[self.mode]

    def generate_scoring_objects(self, bgcs, spectra, tsne_id):
        """
        This method handles constructing the set of objects to be used as input to
        the nplinker "get_links" scoring method, based on the type of objects selected
        from the plots and the currently selected TSNE projection. 

        For most scoring modes it doesn't need to do anything much, but for the
        GCF to Spec and GCF to MolFam modes, the initial selection of BGCs 
        may reference GCFs that don't appear in the current projection 
        (as in these modes we are saying "select ALL GCFs that contain ANY of 
        the selected BGCs"). In this case the unavailable GCFs need to be filtered
        out. 
        """
        scoring_objs = set()

        if self.mode == SCO_MODE_BGC_SPEC or self.mode == SCO_MODE_BGC_MOLFAM:
            # in either of these two modes, we simply want to select the set 
            # of "parent" GCFs for each of the selected BGCs, which is simple:
            scoring_objs = [bgc.parent for bgc in bgcs]
        elif self.mode == SCO_MODE_GCF_SPEC or self.mode == SCO_MODE_GCF_MOLFAM:
            # this pair of modes are more complex. first build a set of *all* 
            # GCFs that contain *any* of the selected BGCs

            # this is the set of GCFs available in the current TSNE projection
            available_gcfs = self.nh.available_gcfs[tsne_id]

            # for each selected BGC...
            for bgc in bgcs:
                if bgc not in self.nh.bgc_gcf_lookup:
                    print('WARNING: missing BGC in lookup: {}'.format(bgc))
                    continue

                # ... get the set of GCFs of which it is a member...
                gcfs_for_bgc = self.nh.bgc_gcf_lookup[bgc]

                # ... and filter out any not in the current TSNE projection
                gcfs_for_bgc = [gcf for gcf in gcfs_for_bgc if gcf in available_gcfs]

                scoring_objs.update(gcfs_for_bgc)
            print('Selection of {} BGCs => {} GCFs'.format(len(bgcs), len(scoring_objs)))
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

def get_radius(datasource):
    # this aims to set a sensible initial radius based on the dimensions of the plot
    minx = min(datasource.data['x'])
    maxx = max(datasource.data['x'])
    return (maxx - minx) / POINT_FACTOR

class ZoomHandler(object):
    """
    This class deals with scaling the circle glyphs on the scatterplots so that
    they stay roughly the same size during zooming. 
    """

    def __init__(self, plot, renderer):
        self.lod_mode = False
        self.renderer = renderer
        self.ratio = POINT_RATIO

        # want to listen for changes to the start/end of the x/y ranges displayed by the plot...
        # can probably get away with only listening for one of these? (and y events seem to be
        # delivered after x)
        # plot.x_range.on_change('end', self.x_cb)
        plot.y_range.on_change('end', self.y_cb)

        # ... and also entering/leaving "level of detail" mode, a lower-detail mode activated
        # during panning/zooming etc. this allows changes producing by panning which we don't
        # care about because it doesn't alter the zoom level
        plot.on_event(LODStart, self.lod_event)
        plot.on_event(LODEnd, self.lod_event)
        plot.on_event(Reset, self.reset_event)

        self.plot = plot

    def reset_event(self, e):
        radius = get_radius(self.renderer.data_source)
        print('reset to  {}'.format(radius))
        self.renderer.glyph.radius = radius

    def update_ratio(self):
        # handle init state when x_/y_range are None
        if self.plot.y_range.start is not None:
            self.renderer.glyph.radius = self.ratio * min(self.plot.x_range.end - self.plot.x_range.start, self.plot.y_range.end - self.plot.y_range.start)

    def x_cb(self, attr, old, new):
        """ 
        Callback for changes to the x range. attr = 'start'/'end', 'old' and 'new'
        give the previous and current values.
        """
        self.update_ratio()

    def y_cb(self, attr, old, new):
        """ 
        Callback for changes to the y range. attr = 'start'/'end', 'old' and 'new'
        give the previous and current values.
        """
        self.update_ratio()

    def lod_event(self, event):
        self.lod_mode = isinstance(event, LODStart)

class CurrentResults(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self.gcfs = []
        self.spec = []
        self.links = {} 
        self.shared_strains = {}

    def update(self, gcfs, spec, links, shared_strains):
        self.gcfs = gcfs
        self.spec = spec
        self.links = links
        self.shared_strains = shared_strains

class NPLinkerBokeh(object):

    def __init__(self, helper):
        self.nh = helper

        # default to TSNE with all bgcs
        self.bgc_tsne_id = 'allbgcs'
        self.bgc_tsne_id_list = list(self.nh.bgc_data.keys())

        self.bgc_datasource = ColumnDataSource(data=self.nh.bgc_data[self.bgc_tsne_id])
        self.spec_datasource = ColumnDataSource(data=self.nh.spec_data)
        # hide singletons by default
        self.singleton_indices = [i for i in range(len(self.spec_datasource.data['x'])) if self.spec_datasource.data['family'][i] != '-1']
        self.spec_datasource_view = CDSView(source=self.spec_datasource, filters=[IndexFilter(self.singleton_indices)])

        self.bgc_zoom = None
        self.spec_zoom = None

        self.score_helper = ScoringHelper(self.nh)

        self.results = CurrentResults()
        self.filter_no_shared_strains = True

    @property
    def bgc_data(self):
        return self.nh.bgc_data[self.bgc_tsne_id]

    @property
    def nplinker(self):
        return self.nh.nplinker

    @property
    def spec_data(self):
        return self.nh.spec_data

    @property
    def spec_indices(self):
        return self.nh.spec_indices

    @property
    def bgc_indices(self):
        return self.nh.bgc_indices[self.bgc_tsne_id]

    def create_plots(self):
        # create the BGC figure and populate it from the datasource (note the dict keys are used in the 
        # call to circle to identify attributes)

        radius = get_radius(self.bgc_datasource)

        f_bgc = figure(tools=TOOLS, toolbar_location='above', title="BGCs (n={})".format(len(self.bgc_data['x'])), sizing_mode='scale_width', name="fig_bgc", output_backend='webgl')
        r_bgc = f_bgc.circle('x', 'y', source=self.bgc_datasource, name=self.bgc_tsne_id,
                radius=radius,
                radius_dimension='max',
                fill_alpha=1.0, 
                fill_color='fill',
                line_color=None,
                selection_color='fill',
                selection_fill_alpha=0.9,
                selection_line_color=None,
                nonselection_fill_alpha=0.0, 
                nonselection_fill_color='#333333', 
                nonselection_line_color=None, 
                nonselection_line_alpha=0)

        self.bgc_zoom = ZoomHandler(f_bgc, r_bgc)

        hover_callback_code = """
            var indices = cb_data.index['1d'].indices;
            if (indices.length > 0)
            {
                $(output_div_id).html('<small>' + multi_prefix + indices.length + multi_suffix + '</small>');
            } 
            else
            {
                $(output_div_id).html('<small>' + empty_text + '</small>');
            }
        """

        callback = CustomJS(args={'output_div_id': '#bgc_hover_info', 
                                    'multi_prefix': 'Click to select ', 
                                    'multi_suffix' : ' BGCs', 
                                    'empty_text' : '(nothing under cursor)'}, 
                                    code=hover_callback_code)
        hover_bgc = HoverTool(tooltips=None, callback=callback)
        f_bgc.add_tools(hover_bgc)

        # create the MolFam figure in the same way
        radius = 0.1
        f_spec = figure(tools=TOOLS, toolbar_location='above', title="Spectra (n={})".format(len(self.nh.spec_data['x'])), sizing_mode='scale_width', name="fig_spec", output_backend='webgl')
        r_spec = f_spec.circle('x', 'y', source=self.spec_datasource, 
                            view=self.spec_datasource_view,
                            fill_alpha=1.0, 
                            radius=radius,
                            radius_dimension='max',
                            fill_color='fill',
                            line_color=None,
                            selection_color='fill',
                            selection_fill_alpha=0.9,
                            selection_line_color=None,
                            nonselection_fill_alpha=0.0, 
                            nonselection_fill_color='#333333', 
                            nonselection_line_color=None, 
                            nonselection_line_alpha=0)

        self.spec_zoom = ZoomHandler(f_spec, r_spec)


        # hover_spec = HoverTool(tooltips=[('index', '$index'), ('Spectra name', '@name')])
        callback = CustomJS(args={'output_div_id': '#spec_hover_info', 
                                    'multi_prefix': 'Click to select ', 
                                    'multi_suffix' : ' spectra', 
                                    'empty_text' : '(nothing under cursor)'}, 
                                    code=hover_callback_code)
        hover_spec = HoverTool(tooltips=None, callback=callback)
        f_spec.add_tools(hover_spec)

        self.fig_bgc = f_bgc
        self.fig_spec = f_spec
        self.ren_bgc = r_bgc
        self.ren_spec = r_spec
        self.hover_bgc = hover_bgc
        self.hover_spec = hover_spec

        self.ds_bgc = self.ren_bgc.data_source
        self.ds_spec = self.ren_spec.data_source

        # self.fig_bgc.on_event('tap', lambda e: self.plot_tapped(self.fig_bgc))
        # self.fig_spec.on_event('tap', lambda e: self.plot_tapped(self.fig_spec))

        self.set_inactive_plot(self.fig_spec)

        return f_bgc, r_bgc, f_spec, r_spec

    def set_tap_behavior(self, plot, btype):
        # this basically prevents clicks from selecting anything on the plot
        taptool = plot.select(type=TapTool)[0]
        taptool.behavior = btype

    def set_inactive_plot(self, inactive):
        print('set inactive')
        if self.fig_bgc == inactive:
            if 'indices' in self.ds_bgc.selected._callbacks and len(self.ds_bgc.selected._callbacks['indices']) > 0:
                self.ds_bgc.selected.remove_on_change('indices', self.bgc_selchanged)
            self.clear_selections()
            self.fig_spec.background_fill_color = ACTIVE_BG_FILL
            self.fig_spec.background_fill_alpha = ACTIVE_BG_ALPHA
            self.fig_bgc.background_fill_color = INACTIVE_BG_FILL
            self.fig_bgc.background_fill_alpha = INACTIVE_BG_ALPHA
            self.ds_spec.selected.on_change('indices', self.spec_selchanged)

            
            self.set_tap_behavior(self.fig_bgc, 'inspect')
            self.set_tap_behavior(self.fig_spec, 'select')

        else:
            if 'indices' in self.ds_spec.selected._callbacks and len(self.ds_spec.selected._callbacks['indices']) > 0:
                self.ds_spec.selected.remove_on_change('indices', self.spec_selchanged)

            self.clear_selections()
            self.fig_bgc.background_fill_color = ACTIVE_BG_FILL
            self.fig_bgc.background_fill_alpha = ACTIVE_BG_ALPHA
            self.fig_spec.background_fill_color = INACTIVE_BG_FILL
            self.fig_spec.background_fill_alpha = INACTIVE_BG_ALPHA
            self.ds_bgc.selected.on_change('indices', self.bgc_selchanged)

            self.set_tap_behavior(self.fig_bgc, 'select')
            self.set_tap_behavior(self.fig_spec, 'inspect')

    # messing with the tools for the plots seems to produce errors in browser console
    # so these aren't currently used
    def enable_bgc_hover(self):
        if self.hover_bgc not in self.fig_bgc.tools:
            self.fig_bgc.add_tools(self.hover_bgc)
            # self.fig_bgc.add_tools(self.taptool_bgc)

    def enable_spec_hover(self):
        if self.hover_spec not in self.fig_spec.tools:
            self.fig_spec.add_tools(self.hover_spec)
            # self.fig_spec.add_tools(self.taptool_spec)

    def disable_bgc_hover(self):
        if self.hover_bgc in self.fig_bgc.tools:
            # self.fig_bgc.tools.active_tap = None
            # self.fig_bgc.tools.active_inspect = None
            self.fig_bgc.tools.remove(self.hover_bgc)
            # self.fig_bgc.tools.remove(self.taptool_bgc)

    def disable_spec_hover(self):
        if self.hover_spec in self.fig_spec.tools:
            # self.fig_spec.tools.active_tap = None
            # self.fig_spec.tools.active_inspect = None
            self.fig_spec.tools.remove(self.hover_spec)
            # self.fig_spec.tools.remove(self.taptool_spec)

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
    TMPL = '''
            <div class="card">
                <div class="card-header" id="{}" style="background-color: #{}">
                    <button class="btn btn-link" type="button" data-toggle="collapse" data-target="#{}">
                        {}
                    </button>
                </div>
                <div id="{}" class="collapse" data-parent="{}">
                    <div class="card-body">
                        {}
                    </div>
                </div>
            </div>'''


    # same as above but with an extra pair of parameters for the "onclick" handler and button id
    TMPL_ON_CLICK = '''
        <div class="card">
            <div class="card-header" id="{}" style="background-color: #{}">
                <button class="btn btn-link" type="button" data-toggle="collapse" data-target="#{}" onclick="{}" id="{}">
                    {}
                </button>
            </div>
            <div id="{}" class="collapse" data-parent="{}">
                <div class="card-body">
                    {}
                </div>
            </div>
        </div>'''

    def gen_gcf_card(self, gcf, index):
        pass

    # TODO
    # this pair of functions could really use some refactoring!

    def update_bgc_output(self):
        """
        Generate the HTML to display a list of BGC objects
        """
        print('>> update_bgc_output for mode {}'.format(self.score_helper.mode_name))
        content = ''

        if len(self.ds_bgc.selected.indices) == 0 and len(self.ds_spec.selected.indices) == 0:
            self.bgc_div.text = content
            return

        # scoring from BGCs to Spectra
        # TODO should split this up to show GCFs
        if self.score_helper.mode == SCO_MODE_BGC_SPEC or self.score_helper.mode == SCO_MODE_GCF_SPEC:
            # in this mode, currently just want to show a simple flat list
            # of the selected BGCs, with the ability to view GCF info as well
            # TODO possibly show grouped by GCF instead? 

            content += '<h4>{} GCFs ({} BGCs)</h4>'.format(len(self.results.gcfs), len(self.ds_bgc.selected.indices))

            # have to construct list of BGCs from the GCFs
            bgcs = set()
            for i, gcf in enumerate(self.results.gcfs):
                bgcs.update(gcf.bgc_list)

            for i, bgc in enumerate(list(bgcs)):
                gcf = bgc.parent

                # handle BGCs not in the current TSNE
                if bgc.name not in self.bgc_indices:
                    print('Skipping BGC {}'.format(bgc))
                    continue

                hdr_id = 'bgc_header_{}'.format(i)
                body_id = 'bgc_body_{}'.format(i)
                title = 'BGC #{}, name={}'.format(i, bgc.name)
                body = '<dl>'
                for attr in ['name', 'strain']:
                    body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(bgc, attr))
                body += '<dt>Parent GCF:</dt><dd>{}</dd>'.format(gcf)
                body += '</dl>'

                content += self.TMPL.format(hdr_id, 'b4c4e8', body_id, title, body_id, 'accordionBGC', body)

        # scoring from Spectra to GCFs
        elif self.score_helper.mode == SCO_MODE_SPEC_GCF:
            uniq_gcfs = set()
            i = 0
            unsorted = []

            if self.results.links is not None:
                # results.links will be a dict keyed by Spectrum objects, values will be
                # lists of (GCF, score) pairs

                # want to construct a 2-level nested list, with the top level entries corresponding to
                # the selected spectra and the contents of each being a list of the GCFs linked to that spectrum

                # for each spectrum selected
                for spec, gcf_scores in self.results.links.items():
                    # construct the object info card header info
                    hdr_id = 'spec_result_header_{}'.format(i)
                    body_id = 'spec_result_body_{}'.format(i)
                    title = '{} GCFs linked to {}'.format(len(gcf_scores), spec)
                    body = '<div class="accordion" id="accordion_spec_{}>'.format(i)
                    # now generate the body content of the card by iterating over the GCFs it matched
                    for j, gcf_score in enumerate(gcf_scores):
                        gcf, score = gcf_score

                        # if this test fails, it means this Spectrum:GCF pair had no shared strains
                        # and will be excluded 
                        if (spec, gcf) not in self.results.shared_strains:
                            continue

                        shared_strains = self.results.shared_strains[(spec, gcf)]
                        uniq_gcfs.add(gcf)

                        # construct a similar object info card header for this GCF
                        gcf_hdr_id = 'gcf_result_header_{}_{}'.format(i, j)
                        gcf_body_id = 'gcf_body_{}_{}'.format(i, j)
                        gcf_title = 'GCF(id={}), score=<strong>{}</strong>, shared strains=<strong>{}</strong>'.format(gcf.id, score, len(shared_strains))
                        gcf_body = '<dl>'
                        for attr in ['id', 'short_gcf_id']:
                            gcf_body += '<dt>{}:</dt> <dd>{}</dd>'.format(attr, getattr(gcf, attr))

                        # add strain information
                        gcf_body += '<dt>shared ({}, {:.2f}%):</dt> <dd>{}</dd>'.format(len(shared_strains), 
                                                                                                    100 * (len(shared_strains) / len(gcf.bgc_list)),
                                                                                                    ', '.join(shared_strains))
                        gcf_body += '<dt>all ({}):</dt> <dd>{}</dd>'.format(len(gcf.bgc_list), ', '.join([bgc.name for bgc in gcf.bgc_list]))
                        gcf_body += '</dl>'
                        body += self.TMPL.format(gcf_hdr_id, 'ffe0b5', gcf_body_id, gcf_title, gcf_body_id, 'accordion_spec_{}'.format(i), gcf_body)
                        
                    body += '</div>'
                    unsorted.append((len(gcf_scores), self.TMPL.format(hdr_id, 'b4c4e8', body_id, title, body_id, 'accordionBGC', body)))
                    i += 1

                content += '<h4>Found {} linked GCFs</h4>'.format(len(uniq_gcfs))
                # sort by number of GCFs linked to each spectra
                sorted_gcfs = sorted(unsorted, key=lambda x: -x[0])
                for _, text in sorted_gcfs:
                    content += text

        self.bgc_div.text = content
        print('<< update_bgc_output')

    def update_spec_output(self):
        """
        Generate the HTML to display a list of Spectrum objects
        """
        print('>> update_spec_output')
        content = ''

        if len(self.ds_bgc.selected.indices) == 0 and len(self.ds_spec.selected.indices) == 0:
            self.spec_div.text = content
            return

        # scoring from Spectra to GCFs
        if self.score_helper.mode == SCO_MODE_SPEC_GCF:
            content += '<h4>{} selected spectra</h4>'.format(len(self.ds_spec.selected.indices))
            i = 0
            # here just need to show a flat list of spectra
            for spec in self.results.spec:
                hdr_id = 'spec_header_{}'.format(i)
                body_id = 'spec_body_{}'.format(i)
                title = 'id={}, family={}'.format(spec.id, spec.family)
                body = '<dl>'
                for attr in ['id', 'spectrum_id', 'family', 'strain_list']:
                    body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(spec, attr))
                # TODO spectrum plot etc here too
                body += '</dl>'
                content += self.TMPL.format(hdr_id, 'adeaad', body_id, title, body_id, 'accordionSpec', body)
                i += 1
            
        # scoring from BGCs to Spectra
        elif self.score_helper.mode == SCO_MODE_BGC_SPEC or self.score_helper.mode == SCO_MODE_GCF_SPEC:
            content += '<h4>{} linked spectra</h4>'.format(len(self.ds_spec.selected.indices))
            i = 0
            unsorted = []

            if self.results.links is not None:
                # results.links will be a dict keyed by Spectrum objects, values will be
                # lists of (GCF, score) pairs

                # want to construct a 2-level nested list, with the top level entries corresponding to
                # the selected spectra and the contents of each being a list of the GCFs linked to that spectrum

                # for each GCF selected
                for gcf, spec_scores in self.results.links.items():
                    # construct the object info card header info
                    hdr_id = 'gcf_result_header_{}'.format(i)
                    body_id = 'gcf_result_body_{}'.format(i)
                    body = '<div class="accordion" id="accordion_gcf_{}">'.format(i)
                    spec_count = 0
                    # now generate the body content of the card by iterating over the spectra it matched
                    for j, spec_score in enumerate(spec_scores):
                        spec, score = spec_score

                        # if this test fails, it means this Spectrum:GCF pair had no shared strains
                        # and will be excluded 
                        if (spec, gcf) not in self.results.shared_strains:
                            continue

                        spec_count += 1
                        shared_strains = self.results.shared_strains[(spec, gcf)]

                        # TODO this could all be neater/templated or something

                        spec_hdr_id = 'spec_result_header_{}_{}'.format(i, j)
                        spec_body_id = 'spec_body_{}_{}'.format(i, j)
                        spec_title = 'Spectrum(id={}), score=<strong>{}</strong>, shared strains=<strong>{}</strong>'.format(spec.id, score, len(shared_strains))
                        spec_body = '<dl>'
                        for attr in ['id', 'spectrum_id', 'family']:
                            spec_body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(spec, attr))
                        
                        # add strain information
                        spec_body += '<dt>shared ({}, {:.2f}%):</dt> <dd>{}</dd>'.format(len(shared_strains), 
                                                                                                    100 * (len(shared_strains) / len(spec.strain_list)),
                                                                                                    ', '.join(shared_strains))
                        spec_body += '<dt>all ({}):</dt> <dd>{}</dd>'.format(len(spec.strain_list), ', '.join(spec.strain_list))
                        spec_body += '</dl>'
                        
                        # set up the chemdoodle plot so it appears when the entry is expanded 
                        spec_btn_id = 'spec_btn_{}_{}'.format(i, j)
                        spec_plot_id = 'spec_plot_{}_{}'.format(i, j)
                        spec_body += '<canvas id="{}"></canvas>'.format(spec_plot_id)
                        # note annoying escaping required here, TODO better way of doing this?
                        spec_onclick = 'setupPlot(\'{}\', \'{}\', \'{}\');'.format(spec_btn_id, spec_plot_id, spec.to_jcamp_str())

                        body += self.TMPL_ON_CLICK.format(spec_hdr_id, 'ffe0b5', spec_body_id, spec_onclick, spec_btn_id, spec_title, spec_body_id, 'accordion_gcf_{}'.format(i), spec_body)

                    body += '</div>'
                    title = '{} spectra linked to GCF(id={})'.format(spec_count, gcf.id)
                    unsorted.append((len(spec_scores), self.TMPL.format(hdr_id, 'adeaad', body_id, title, body_id, 'accordionSpec', body)))
                    i += 1

            for _, text in sorted(unsorted, key=lambda x: -x[0]):
                content += text

        self.spec_div.text = content
        print('<< update_spec_output')

    def get_links(self):
        """
        Handles the scoring process
        """
        nplinker = self.nh.nplinker
        sel_indices = []

        # first check the scoring mode, which is determined by the state 
        # of the two sets of UI controls. based on the current mode, 
        # extract the appropriate set of selected indices from the plot
        bgcs = []
        spectra = []
        if self.score_helper.from_genomics:
            sel_indices = self.ds_bgc.selected.indices
            for i in range(len(sel_indices)):
                bgcs.append(nplinker.lookup_bgc(self.bgc_data['name'][sel_indices[i]]))
        else:
            sel_indices = self.ds_spec.selected.indices
            for i in range(len(sel_indices)):
                spectra.append(nplinker.lookup_spectrum(self.spec_data['name'][sel_indices[i]]))

        mode = self.score_helper.mode

        if len(sel_indices) == 0:
            self.update_alert('No objects selected', 'danger')
            return

        # TODO why did i do this
        nplinker.clear_links() # TODO hack 
        self.results.clear()
        current_method = nplinker.scoring.enabled()[self.scoring_method_group.active]

        # obtain a list of objects to be passed to the scoring function
        scoring_objs = self.score_helper.generate_scoring_objects(bgcs, spectra, self.bgc_tsne_id)
        if scoring_objs is None:
            self.update_alert('ERROR: something went wrong in scoring!', 'danger')
            self.clear_selections()
            return

        # this method returns the input objects that have links based on the selected scoring parameters
        # (ie it will be a (possibly complete, possibly empty) subset of the original list
        objs_with_links = nplinker.get_links(scoring_objs, scoring_method=current_method)
        if len(objs_with_links) > 0:

            # now we know which objects actually have links, want to obtain the 
            # objects at the other end of those links

            # the type parameter limits the links returned by links_for_obj() to objects
            # of the supplied class. Otherwise for example when scoring GCFs, it will
            # return both Spectra and MolFams, and don't want to do that here
            t = Spectrum
            if mode == SCO_MODE_BGC_MOLFAM or mode == SCO_MODE_GCF_MOLFAM:
                t = MolecularFamily
            elif mode == SCO_MODE_SPEC_GCF or mode == SCO_MODE_MOLFAM_GCF:
                t = GCF

            # create a dict with keys = input objects, values = link info for that object
            objs_with_scores = {}
            for link_obj in objs_with_links:
                objs_with_scores[link_obj] = nplinker.links_for_obj(link_obj, current_method, type_=t)

            # if scoring from Spectra/MolFam to GCF, now need to walk through the list of 
            # results and exclude any which feature GCFs that do not appear in the current
            # TSNE projection (this doesn't need done for the SCO_MODE_BGC_ and SCO_MODE_GCF_
            # modes because in the former case only available GCFs will have been selected in
            # the first place, and in the latter case unavailable GCFs have been filtered
            # out already by calling generate_scoring_objects())
            if mode == SCO_MODE_SPEC_GCF or mode == SCO_MODE_MOLFAM_GCF:
                print('Excluding GCFs from objs_with_scores')
                for spec_or_fam in list(objs_with_scores.keys()):
                    gcfs = [(gcf, score) for gcf, score in objs_with_scores[spec_or_fam] if gcf in self.nh.available_gcfs[self.bgc_tsne_id]]

                    # if this result now happens to have zero GCFs, remove it completely
                    if len(gcfs) == 0:
                        del objs_with_scores[spec_or_fam]
                        objs_with_links.remove(spec_or_fam)
                        print('Removing result {}'.format(spec_or_fam))
                    else:
                        # otherwise just overwrite the original list
                        objs_with_scores[spec_or_fam] = gcfs

            # now want to retrieve shared strain information for all these objects. the filter_no_shared_strains 
            # parameter determines whether results will be returned for links with no shared strains or not. 
            otherobjs = set()
            for link_obj, sco_data in objs_with_scores.items():
                otherobjs.update(sco_obj for sco_obj, sco in sco_data)
            all_shared_strains = self.nh.nplinker.get_common_strains(objs_with_links, list(otherobjs), self.filter_no_shared_strains)

            # the all_shared_strains dict is indexed by object pairs, ALWAYS with Spectra/MolFam
            # as the first entry in the tuple. the next step in filtering the results 
            # (if self.filter_no_shared_strains is True) is to go through the results and remove links where
            # that pair does not appear in all_shared_strains
            if self.filter_no_shared_strains:
                print('Filtering out results with no shared strains, initial count={}'.format(len(objs_with_scores[objs_with_links[0]])))
                print(objs_with_links)

                # for each object which has links
                # TODO can these be combined in a better way?
                if mode == SCO_MODE_BGC_SPEC or mode == SCO_MODE_GCF_SPEC or mode == SCO_MODE_BGC_MOLFAM or mode == SCO_MODE_GCF_MOLFAM:
                    for link_obj, sco_objs_and_scores in objs_with_scores.items():
                        to_keep = set()
                        for sco_obj, score in sco_objs_and_scores:
                            # check if each pairing has any shared strains and if so add to set
                            if (sco_obj, link_obj) in all_shared_strains:
                                to_keep.add(sco_obj)

                        # rebuild list including only those in the set
                        objs_with_scores[link_obj] = [(sco_obj, score) for sco_obj, score in sco_objs_and_scores if sco_obj in to_keep]
                else:
                    for link_obj, sco_objs_and_scores in objs_with_scores.items():
                        to_keep = set()
                        for sco_obj, score in sco_objs_and_scores:
                            # check if each pairing has any shared strains and if so add to set
                            if (link_obj, sco_obj) in all_shared_strains:
                                to_keep.add(sco_obj)

                        # rebuild list including only those in the set
                        objs_with_scores[link_obj] = [(sco_obj, score) for sco_obj, score in sco_objs_and_scores if sco_obj in to_keep]

                print('After filtering out results with no shared strains, initial count={}'.format(len(objs_with_scores[objs_with_links[0]])))

            # final step here is to take the remaining list of scoring objects and map them
            # back to indices on the plot so they can be highlighted
            selected = set()
            # TODO is this really needed???
            gcfs = set()
            for link_obj in objs_with_links:
                score_obj_indices = []
                if mode == SCO_MODE_BGC_SPEC or mode == SCO_MODE_GCF_SPEC:
                    # for these modes, we can just lookup spectrum indices directly
                    score_obj_indices = [self.spec_indices[spec.spectrum_id] for (spec, score) in objs_with_scores[link_obj]]
                    gcfs.add(link_obj)
                elif mode == SCO_MODE_SPEC_GCF or mode == SCO_MODE_MOLFAM_GCF:
                    # here we have a list of GCFs, but the plot shows BGCs. so in order
                    # to highlight the BGCs, need to iterate over the list of GCFs 
                    score_obj_indices = set()
                    for (gcf, score) in objs_with_scores[link_obj]:
                        bgc_names = [bgc.name for bgc in gcf.bgc_list]
                        for n in bgc_names:
                            try:
                                score_obj_indices.add(self.bgc_indices[n])
                            except KeyError:
                                print('Warning: missing index for BGC: {} in GCF: {}'.format(n, gcf))
                    score_obj_indices = list(score_obj_indices)

                # add these indices to the overall set
                selected.update(score_obj_indices)

        print('After filtering, remaining objs={}'.format(len(objs_with_links)))

        if len(objs_with_links) > 0:
            # take the indexes of the corresponding spectra from the datasource and highlight on the other plot
            if mode == SCO_MODE_BGC_SPEC or mode == SCO_MODE_GCF_SPEC:
                self.ds_spec.selected.indices = list(selected)
                self.results.update(gcfs, None, objs_with_scores, all_shared_strains)
            elif mode == SCO_MODE_SPEC_GCF:
                self.ds_bgc.selected.indices = list(selected)
                self.results.update(None, spectra, objs_with_scores, all_shared_strains)

            self.update_alert('{} objects with links found'.format(len(objs_with_scores), 'primary'))
        else:
            print('No links found!')
            if mode == SCO_MODE_BGC_SPEC or mode == SCO_MODE_GCF_SPEC:
                self.ds_spec.selected.indices = []
            else:
                self.ds_bgc.selected.indices = []
            self.results.clear()
            self.update_alert('No links found for last selection', 'danger')

        self.update_bgc_output()
        self.update_spec_output()

    def bgc_selchanged(self, attr, old, new):
        """
        Callback for changes to the BGC datasource selected indices
        """
        # if the new selection contains any valid indices (it may be empty)
        if len(new) > 0:
            # get links for the selected BGCs and update the plots
            self.get_links()
        else:
            # if selection is now empty, clear the selection on the spectra plot too
            self.ds_spec.selected.indices = []

            self.update_bgc_output()
            self.update_spec_output()

    def spec_selchanged(self, attr, old, new):
        """
        Callback for changes to the Spectra datasource selected indices
        """
        # if the new selection contains any valid indices (it may be empty)
        if len(new) > 0:
            # get links for the selected spectra and update the plots
            self.get_links()
        else:
            # if selection is now empty, clear the selection on the BGC plot too
            self.ds_bgc.selected.indices = []

            self.update_bgc_output()
            self.update_spec_output()

    def metcalf_percentile_callback(self, attr, old, new):
        self.nh.nplinker.scoring.metcalf.sig_percentile = new
        self.get_links()

    def likescore_cutoff_callback(self, attr, old, new):
        self.nh.nplinker.scoring.likescore.cutoff = new / 100.0
        self.get_links()

    def hg_prob_callback(self, attr, old, new):
        self.nh.nplinker.scoring.hg.prob = new / 100.0
        self.get_links()

    def scoring_method_callback(self, attr, old, new):
        self.get_links()

    def clear_selections(self):
        self.ds_bgc.selected.indices = []
        self.ds_spec.selected.indices = []

    def get_scoring_mode_text(self):
        return SCO_MODE_NAMES[self.score_helper.mode]

    def sco_mode_changed(self):
        self.update_alert('Scoring mode is now <strong>{}</strong>'.format(self.get_scoring_mode_text()))

    # NOT CURRENTLY USED
    # def plot_tapped(self, plot):
    #     if plot == self.fig_bgc:
    #         print('TAPPED BGC')
    #         self.score_helper.set_genomics()
    #         self.set_inactive_plot(self.fig_spec)
    #     else:
    #         print('TAPPED SPEC')
    #         self.score_helper.set_metabolomics()
    #         self.set_inactive_plot(self.fig_bgc)
    #     self.sco_mode_changed()

    def mb_mode_callback(self, attr, old, new):
        self.score_helper.update_metabolomics(new)
        self.set_inactive_plot(self.fig_bgc)
        self.sco_mode_changed()
        self.update_plot_select_state(False)
    
    def ge_mode_callback(self, attr, old, new):
        self.score_helper.update_genomics(new)
        self.set_inactive_plot(self.fig_spec)
        self.sco_mode_changed()
        self.update_plot_select_state(True)

    def plot_toggles_callback(self, attr, old, new):
        self.ren_bgc.glyph.fill_alpha = 0.6 if PLOT_ALPHA in new else 1.0
        self.ren_spec.glyph.fill_alpha = 0.6 if PLOT_ALPHA in new else 1.0

        self.ren_bgc.glyph.fill_color = 'fill' if PLOT_CMAP in new else '#449944'
        self.ren_spec.glyph.fill_color = 'fill' if PLOT_CMAP in new else '#444499'

        if PLOT_SINGLETONS in new:
            self.ren_spec.view.filters = [IndexFilter(list(range(len(self.spec_datasource.data['x']))))]
        else:
            self.ren_spec.view.filters = [IndexFilter(self.singleton_indices)]

        if PLOT_PRESERVE_COLOUR in new:
            self.ren_bgc.selection_glyph.fill_color = 'fill'
            self.ren_spec.selection_glyph.fill_color = 'fill'
        else:
            self.ren_bgc.selection_glyph.fill_color = '#ff7f00'
            self.ren_spec.selection_glyph.fill_color = '#ff7f00'

        self.filter_no_shared_strains = PLOT_ONLY_SHARED_STRAINS in new

    def tsne_id_callback(self, attr, old, new):
        self.clear_selections()
        self.bgc_tsne_id = new
        self.update_alert('TSNE set to {}'.format(new), 'primary')
        # messy but it works and seems to be best way of doing this
        newdata = self.nh.bgc_data[new]
        self.bgc_datasource.data.update(
                                            x=newdata['x'],
                                            y=newdata['y'],
                                            name=newdata['name'],
                                            strain=newdata['strain'],
                                            gcf=newdata['gcf'],
                                            fill=newdata['fill'])

    def update_alert(self, msg, alert_class='primary'):
        self.alert_div.text = '<div class="alert alert-{}" role="alert">{}</div>'.format(alert_class, msg)

    def clear_alert(self):
        self.alert_div.text = ''

    def update_plot_select_state(self, gen_to_met):
        if gen_to_met:
            self.plot_select.label = SCORING_GEN_TO_MET
            self.plot_select.css_classes = ['button-genomics']
        else:
            self.plot_select.label = SCORING_MET_TO_GEN
            self.plot_select.css_classes = ['button-metabolomics']

    def plot_select_callback(self, val):
        if val:
            self.score_helper.set_genomics()
            self.set_inactive_plot(self.fig_spec)
        else:
            self.score_helper.set_metabolomics()
            self.set_inactive_plot(self.fig_bgc)
        self.update_plot_select_state(val)

    def bokeh_layout(self):
        self.spec_div = Div(text="", sizing_mode='scale_height', name='spec_div')
        self.bgc_div = Div(text="", sizing_mode='scale_height', name='bgc_div')

        self.ds_bgc.selected.on_change('indices', self.bgc_selchanged)

        self.plot_toggles = CheckboxGroup(active=[PLOT_CMAP, PLOT_PRESERVE_COLOUR, PLOT_ONLY_SHARED_STRAINS], labels=PLOT_TOGGLES, name='plot_toggles')
        self.plot_toggles.on_change('active', self.plot_toggles_callback)

        self.tsne_id_select = Select(title='BGC TSNE:', value=self.bgc_tsne_id, options=self.bgc_tsne_id_list, name='tsne_id_select')
        self.tsne_id_select.on_change('value', self.tsne_id_callback)

        self.mb_mode_group = RadioGroup(labels=METABOLOMICS_SCORING_MODES, active=0, name='mb_mode_group', css_classes=['modes-metabolomics'])
        self.mb_mode_group.on_change('active', self.mb_mode_callback)

        self.plot_select = Toggle(label=SCORING_GEN_TO_MET, name='plot_select', active=True, css_classes=['button-genomics'])
        self.plot_select.on_click(self.plot_select_callback)

        self.ge_mode_group = RadioGroup(labels=GENOMICS_SCORING_MODES, active=0, name='ge_mode_group')
        self.ge_mode_group.on_change('active', self.ge_mode_callback)

        self.scoring_method_group = RadioGroup(labels=[m for m in self.nh.nplinker.scoring.enabled_names()], active=0, name='scoring_method_group')
        self.scoring_method_group.on_change('active', self.scoring_method_callback)

        # metcalf stuff
        self.metcalf_percentile = Slider(start=70, end=100, value=self.nh.nplinker.scoring.metcalf.sig_percentile, step=1, title='[metcalf] sig_percentile')
        self.metcalf_percentile.on_change('value', self.metcalf_percentile_callback)
        # likescore
        self.likescore_cutoff = Slider(start=0, end=100, value=int(100 * self.nh.nplinker.scoring.likescore.cutoff), step=1, title='[likescore] cutoff x 100 = ')
        self.likescore_cutoff.on_change('value', self.likescore_cutoff_callback)
        # hg
        self.hg_prob = Slider(start=0, end=100, value=int(100 * self.nh.nplinker.scoring.hg.prob), step=1, title='[hg] prob x 100 = ')
        self.hg_prob.on_change('value', self.hg_prob_callback)

        active_sliders = []
        if self.nh.nplinker.scoring.metcalf.enabled:
            active_sliders.append(self.metcalf_percentile)
        if self.nh.nplinker.scoring.hg.enabled:
            active_sliders.append(self.hg_prob)
        if self.nh.nplinker.scoring.likescore.enabled:
            active_sliders.append(self.likescore_cutoff)
        self.sliders = row(active_sliders, name='sliders')

        # for debug output etc 
        self.debug_div = Div(text="", name='debug_div')

        # status updates/errors
        self.alert_div = Div(text="", name='alert_div')

        # "reset everything" button
        self.reset_button = Button(name='reset_button', label='Reset state')
        # no python method to reset plots, for some reason...
        self.reset_button.js_on_click(CustomJS(args=dict(fig_bgc=self.fig_bgc, fig_spec=self.fig_spec), 
                                               code=""" 
                                                    fig_bgc.reset.emit();
                                                    fig_spec.reset.emit();
                                                """))

        curdoc().add_root(self.plot_toggles)
        curdoc().add_root(self.tsne_id_select)
        curdoc().add_root(self.reset_button)
        curdoc().add_root(self.alert_div)
        curdoc().add_root(self.fig_spec)
        curdoc().add_root(self.fig_bgc)
        curdoc().add_root(self.spec_div)
        curdoc().add_root(self.bgc_div)
        curdoc().add_root(self.plot_select)
        curdoc().add_root(self.mb_mode_group)
        curdoc().add_root(self.ge_mode_group)
        curdoc().add_root(self.scoring_method_group)
        curdoc().add_root(self.sliders)
        curdoc().add_root(self.debug_div)
        curdoc().title = 'nplinker webapp'

# server_lifecycle.py adds a .nh attr to the current Document instance, use that
# to access the already-created NPLinker instance plus the TSNE data
nb = NPLinkerBokeh(curdoc().nh)
nb.create_plots()
nb.bokeh_layout()
nb.update_alert('Initialised OK!', 'success')

