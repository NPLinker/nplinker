from bokeh.models.widgets import RadioGroup, Slider, Div, Dropdown, CheckboxButtonGroup, Select
from bokeh.layouts import row, column, widgetbox
from bokeh.models import CustomJS
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, HoverTool, TapTool, CDSView, IndexFilter
from bokeh.events import LODEnd, LODStart, Reset

from nplinker import Spectrum, GCF, MolecularFamily

# TOOLS="crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
TOOLS = "crosshair,pan,wheel_zoom,box_zoom,tap,reset,save,box_select"

PW, PH = 500, 500

POINT_RATIO = 1 / 200.0

SCORING_MODES = ['BGCs to Spectra', 'Spectra to GCFs']
SCO_MODE_BGC_SPEC, SCO_MODE_SPEC_BGC = range(2)
SCORING_MODES_ENUM = [SCO_MODE_BGC_SPEC, SCO_MODE_SPEC_BGC]
PLOT_TOGGLES = ['Toggle alpha-blending', 'Toggle colormaps', 'Toggle singletons']
PLOT_ALPHA, PLOT_CMAP, PLOT_SINGLETONS = range(3)
PLOT_TOGGLES_ENUM = [PLOT_ALPHA, PLOT_CMAP, PLOT_SINGLETONS]

def get_radius(datasource):
    # this aims to set a sensible initial radius based on the dimensions of the plot
    minx = min(datasource.data['x'])
    maxx = max(datasource.data['x'])
    return (maxx - minx) / 180.0

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
        self.mode = SCO_MODE_BGC_SPEC

    def clear(self):
        self.bgcs = []
        self.spec = []
        self.links = {} 
        self.shared_strains = {}

    def update(self, bgcs, spec, links, shared_strains):
        self.bgcs = bgcs
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
        self.spec_datasource_view = CDSView(source=self.spec_datasource, filters=[IndexFilter(list(range(len(self.spec_datasource.data['x']))))])

        self.bgc_zoom = None
        self.spec_zoom = None

        self.results = CurrentResults()

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
        return self.nh.bgc_indices

    def create_plots(self):
        # create the BGC figure and populate it from the datasource (note the dict keys are used in the 
        # call to circle to identify attributes)

        radius = get_radius(self.bgc_datasource)

        f_bgc = figure(tools=TOOLS, toolbar_location='above', title="BGCs (n={})".format(len(self.bgc_data['x'])), sizing_mode='scale_width', name="fig_bgc")
        r_bgc = f_bgc.circle('x', 'y', source=self.bgc_datasource, name=self.bgc_tsne_id,
                radius=radius,
                radius_dimension='max',
                fill_alpha=1.0, 
                fill_color='fill',
                line_color=None,
                # selection_color='#ff0000')
                selection_color='#ff7f00', 
                selection_fill_alpha=0.9,
                nonselection_fill_alpha=0.0, 
                nonselection_fill_color='#333333', 
                nonselection_line_color=None, 
                nonselection_line_alpha=0)

        self.bgc_zoom = ZoomHandler(f_bgc, r_bgc)

        hover_callback_code = """
            var indices = cb_data.index['1d'].indices;
            if (indices.length > 0)
            {
                $(output_div_id).text(multi_prefix + indices.length + multi_suffix);
            } 
            else
            {
                $(output_div_id).text(empty_text);
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
        f_spec = figure(tools=TOOLS, toolbar_location='above', title="Spectra (n={})".format(len(self.nh.spec_data['x'])), sizing_mode='scale_width', name="fig_spec")
        r_spec = f_spec.circle('x', 'y', source=self.spec_datasource, 
                            view=self.spec_datasource_view,
                            fill_alpha=1.0, 
                            radius=radius,
                            radius_dimension='max',
                            #fill_color='#449944', 
                            fill_color='fill',
                            line_color=None,
                            selection_color='#ff7f00', 
                            selection_fill_alpha=0.9,
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

        self.fig_bgc.on_event('tap', lambda e: self.scoring_mode_callback_tap(SCO_MODE_BGC_SPEC))
        self.fig_spec.on_event('tap', lambda e: self.scoring_mode_callback_tap(SCO_MODE_SPEC_BGC))

        return f_bgc, r_bgc, f_spec, r_spec

    def enable_bgc_hover(self):
        if self.hover_bgc not in self.fig_bgc.tools:
            self.fig_bgc.add_tools(self.hover_bgc)
            self.fig_bgc.add_tools(self.taptool_bgc)
    
    def enable_spec_hover(self):
        if self.hover_spec not in self.fig_spec.tools:
            self.fig_spec.add_tools(self.hover_spec)
            self.fig_spec.add_tools(self.taptool_spec)
    
    def disable_bgc_hover(self):
        print(self.fig_bgc.tools)
        if self.hover_bgc in self.fig_bgc.tools:
            self.fig_bgc.tools.active_tap = None
            self.fig_bgc.tools.active_inspect = None
            self.fig_bgc.tools.remove(self.hover_bgc)
            self.fig_bgc.tools.remove(self.taptool_bgc)
        print(self.fig_bgc.tools)
    
    def disable_spec_hover(self):
        if self.hover_spec in self.fig_spec.tools:
            self.fig_spec.tools.active_tap = None
            self.fig_spec.tools.active_inspect = None
            self.fig_spec.tools.remove(self.hover_spec)
            self.fig_spec.tools.remove(self.taptool_spec)

    def update_debug_output(self):
        self.debug_div.text = '<pre>{}</pre>'.format(str(self.nh.nplinker.scoring))

    # format params:
    # - id for header, e.g spec_heading_1
    # - colour for header
    # - id for body, e.g. spec_body_1
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

    def update_bgc_output(self):
        """
        Generate the HTML to display a list of BGC objects
        """
        print('>> update_bgc_output for mode {}'.format(self.results.mode))
        if self.results.mode == SCO_MODE_BGC_SPEC:
            # in BGC to Spec mode, just list the BGCs directly
            test = '<h4>{} selected BGCs</h4>'.format(len(self.ds_bgc.selected.indices))
            for i, bgc in enumerate(self.results.bgcs):
                gcf = bgc.parent
                hdr_id = 'bgc_header_{}'.format(i)
                body_id = 'bgc_body_{}'.format(i)
                title = 'BGC #{}, name={}'.format(i, self.bgc_data['name'][i])
                body = '<dl>'
                for attr in ['name', 'strain']:
                    body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(bgc, attr))
                body += '<dt>Parent GCF:</dt><dd>{}</dd>'.format(gcf)
                body += '</dl>'

                test += self.TMPL.format(hdr_id, 'b4c4e8', body_id, title, body_id, 'accordionBGC', body)
        elif self.results.mode == SCO_MODE_SPEC_BGC:
            # in Spec to BGC mode, have to retrieve the BGCs from the .links attr instead
            # and display with scores
            print('Outputting BGC links')
            uniq_gcfs = set()
            test = ''
            if self.results.links is not None:
                # links here will be a dict keyed by Spectrum objects, values will be
                # lists of (GCF, score) pairs, so need to extract the BGCs for each one
                i = 0
                for spec, gcf_scores in self.results.links.items():
                    for j, gcf_score in enumerate(gcf_scores):
                        gcf, score = gcf_score
                        uniq_gcfs.add(gcf)

                        for bgc in gcf.bgc_list:
                            hdr_id = 'bgc_header_{}'.format(i)
                            body_id = 'bgc_body_{}'.format(i)
                            title = 'BGC #{}, score={}, name={}'.format(i, score, bgc.name)
                            body = '<dl>'
                            for attr in ['name', 'strain']:
                                body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(bgc, attr))
                            body += '<dt>Parent GCF:</dt><dd>{}</dd>'.format(gcf)
                            body += '<dt>Score:</dt><dd>{}</dd>'.format(score)
                            body += '</dl>'
                            test += self.TMPL.format(hdr_id, 'b4c4e8', body_id, title, body_id, 'accordionBGC', body)
                            i += 1

            test = '<h4>{} linked BGCs ({} GCFs)</h4>'.format(len(self.ds_bgc.selected.indices), len(uniq_gcfs)) + test

        self.bgc_div.text = test
        print('<< update_bgc_output')

    def update_spec_output(self):
        """
        Generate the HTML to display a list of Spectrum objects
        """
        print('>> update_spec_output')
        i = 0
        if self.results.mode == SCO_MODE_SPEC_BGC:
            test = '<h4>{} selected spectra</h4>'.format(len(self.ds_spec.selected.indices))
            for spec in self.results.spec:
                # TODO
                hdr_id = 'spec_header_{}'.format(i)
                body_id = 'spec_body_{}'.format(i)
                title = 'id={}, family={}'.format(spec.id, spec.family)
                body = '<dl>'
                for attr in ['id', 'spectrum_id', 'family', 'strain_list']:
                    body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(spec, attr))
                # if self.results.shared_strains is not None: # TODO!
                #     body += '<dt>shared strains</dt><dd>{} (with {})</dd>'.format(self.results.shared_strains[spec], gcf)
                body += '</dl>'
                test += self.TMPL.format(hdr_id, 'adeaad', body_id, title, body_id, 'accordionSpec', body)
                i += 1
        elif self.results.mode == SCO_MODE_BGC_SPEC:
            test = '<h4>{} linked spectra</h4>'.format(len(self.ds_spec.selected.indices))
            seen_gcfs = set()
            if self.results.links is not None:
                for gcf, spec_scores in self.results.links.items():
                    for j, spec_score in enumerate(spec_scores):
                        spec, score = spec_score
                        hdr_id = 'spec_header_{}'.format(i)
                        body_id = 'spec_body_{}'.format(i)
                        title = 'id={}, score={}'.format(spec.id, score)
                        body = '<dl>'
                        for attr in ['id', 'spectrum_id', 'family', 'strain_list']:
                            body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(spec, attr))
                        if self.results.shared_strains is not None:
                            for bgc in self.results.bgcs:
                                if (spec, bgc.parent) in self.results.shared_strains and bgc.parent not in seen_gcfs:
                                    if len(self.results.shared_strains[(spec, bgc.parent)]) > 0:
                                        body += '<dt>shared strains</dt><dd>{} (with GCF {}, BGC {})</dd>'.format(self.results.shared_strains[(spec, bgc.parent)], bgc.parent, bgc)
                                    seen_gcfs.add(bgc.parent)
                        body += '</dl>'
                        test += self.TMPL.format(hdr_id, 'adeaad', body_id, title, body_id, 'accordionSpec', body)
                        i += 1

        self.spec_div.text = test
        print('<< update_spec_output')

    def get_links(self):
        """
        Retrieves links for current selected objects, if any.
        """

        # first check the scoring mode
        sel_indices = []
        mode = self.results.mode
        if mode == SCO_MODE_BGC_SPEC:
            sel_indices = self.ds_bgc.selected.indices
        elif mode == SCO_MODE_SPEC_BGC:
            sel_indices = self.ds_spec.selected.indices
        else:
            return

        if len(sel_indices) == 0:
            print('No indices selected')
            return

        print('Selected {} objects'.format(len(sel_indices)))
        nplinker = self.nh.nplinker
        nplinker.clear_links() # TODO hack
        self.results.clear()
        current_method = nplinker.scoring.enabled()[self.scoring_method_group.active]
        print('Scoring method is: {}'.format(current_method))

        scoring_objs = set()
        if mode == SCO_MODE_BGC_SPEC: # indices are into BGC plot
            bgcs = []
            # for each BGC, get its parent GCF
            for i in range(len(sel_indices)):
                bgc = nplinker.lookup_bgc(self.bgc_data['name'][sel_indices[i]])
                scoring_objs.add(bgc.parent)
                bgcs.append(bgc)
        elif mode == SCO_MODE_SPEC_BGC:
            specs = []
            for i in range(len(sel_indices)):
                spec = nplinker.lookup_spectrum(self.spec_data['name'][sel_indices[i]])
                scoring_objs.add(spec)
                specs.append(spec)
        else:
            # TODO other modes
            pass

        scoring_objs = list(scoring_objs)

        print('Initial setup done')
        # this method returns the input objects that have links based on the selected scoring parameters
        objs_with_links = nplinker.get_links(scoring_objs, scoring_method=current_method)
        if len(objs_with_links) > 0:
            objs_with_scores = {}
            selected = set()
            t = Spectrum # default for BGC to Spec
            if mode == SCO_MODE_SPEC_BGC:
                t = GCF
            # TODO other modes

            otherobjs = set()
            for link_obj in objs_with_links:
                objs_with_scores[link_obj] = nplinker.links_for_obj(link_obj, current_method, type_=t)
                if mode == SCO_MODE_BGC_SPEC:
                    score_obj_indices = [self.spec_indices[spec.spectrum_id] for (spec, score) in objs_with_scores[link_obj]]
                elif mode == SCO_MODE_SPEC_BGC:
                    score_obj_indices = set()
                    for (gcf, score) in objs_with_scores[link_obj]:
                        bgc_names = [bgc.name for bgc in gcf.bgc_list]
                        for n in bgc_names:
                            try:
                                score_obj_indices.add(self.bgc_indices[n])
                                # print('Found bgc {}'.format(n))
                            except KeyError:
                                print('Missing index for BGC: {}'.format(n))
                    score_obj_indices = list(score_obj_indices)

                selected.update(score_obj_indices)

                for obj in objs_with_scores[link_obj]:
                    otherobjs.add(obj[0])

                # TODO shared strains
                # shared_strains = nplinker.get_common_strains([gcf], [spec[0] for spec in spectra_and_scores])

            all_shared_strains = self.nh.nplinker.get_common_strains(objs_with_links, list(otherobjs))

            # take the indexes of the corresponding spectra from the datasource and highlight on the other plot
            if mode == SCO_MODE_BGC_SPEC:
                self.ds_spec.selected.indices = list(selected)
                self.results.update(bgcs, None, objs_with_scores, all_shared_strains)
            elif mode == SCO_MODE_SPEC_BGC:
                self.ds_bgc.selected.indices = list(selected)
                self.results.update(None, specs, objs_with_scores, all_shared_strains)

        else:
            print('No links found!')
            if mode == SCO_MODE_BGC_SPEC:
                self.ds_spec.selected.indices = []
            else:
                self.ds_bgc.selected.indices = []
            self.results.clear()
        
        self.update_bgc_output()
        self.update_spec_output()
        self.update_debug_output()

    def bgc_selchanged(self, attr, old, new):
        """
        Callback for changes to the BGC datasource selected indices
        """
        print('bgc_selchanged')
        # if the new selection contains any valid indices (it may be empty)
        if len(new) > 0:
            # get links for the selected BGCs and update the plots
            self.get_links()
        else:
            # if selection is now empty, clear the selection on the spectra plot too
            self.ds_spec.selected.indices = []

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

    def metcalf_percentile_callback(self, attr, old, new):
        self.nh.nplinker.scoring.metcalf.sig_percentile = new
        self.get_links()
        self.update_debug_output()

    def likescore_cutoff_callback(self, attr, old, new):
        self.nh.nplinker.scoring.likescore.cutoff = new / 100.0
        self.get_links()
        self.update_debug_output()

    def hg_prob_callback(self, attr, old, new):
        self.nh.nplinker.scoring.hg.prob = new / 100.0
        self.get_links()
        self.update_debug_output()

    def scoring_method_callback(self, attr, old, new):
        self.get_links()
        self.update_debug_output()

    def scoring_mode_callback_tap(self, newmode):
        self.scoring_mode_group.active = newmode

    def scoring_mode_callback(self, attr, old, new):
        self.results.mode = SCORING_MODES_ENUM[new]
        print('Score mode: {}'.format(self.results.mode))
        if self.results.mode == SCO_MODE_BGC_SPEC:
            # self.enable_bgc_hover()
            # self.disable_spec_hover()
            # only want to listen for selection changes on bgc plot
            self.ds_bgc.selected.on_change('indices', self.bgc_selchanged)
            if len(self.ds_spec.selected._callbacks['indices']) > 0:
                self.ds_spec.selected.remove_on_change('indices', self.spec_selchanged)
        elif self.results.mode == SCO_MODE_SPEC_BGC:
            # self.disable_bgc_hover()
            # self.enable_spec_hover()
            # only want to listen for selection changes on spec plot
            if len(self.ds_bgc.selected._callbacks['indices']) > 0:
                self.ds_bgc.selected.remove_on_change('indices', self.bgc_selchanged)
            self.ds_spec.selected.on_change('indices', self.spec_selchanged)

        # clear existing selections on mode changes
        self.ds_bgc.selected.indices = []
        self.ds_spec.selected.indices = []

    def plot_toggles_callback(self, attr, old, new):
        self.ren_bgc.glyph.fill_alpha = 0.6 if PLOT_ALPHA in new else 1.0
        self.ren_spec.glyph.fill_alpha = 0.6 if PLOT_ALPHA in new else 1.0

        self.ren_bgc.glyph.fill_color = 'fill' if PLOT_CMAP in new else '#449944'
        self.ren_spec.glyph.fill_color = 'fill' if PLOT_CMAP in new else '#444499'

        if PLOT_SINGLETONS in new:
            self.ren_spec.view.filters = [IndexFilter(list(range(len(self.spec_datasource.data['x']))))]
        else:
            indices = [i for i in range(len(self.spec_datasource.data['x'])) if self.spec_datasource.data['family'][i] != '-1']
            # for i in range(len(self.spec_datasource.data['x'])):
            #     if self.spec_datasource.data['family'][i] != '-1':
            #         indices.append(i)

            self.ren_spec.view.filters = [IndexFilter(indices)]

    def tsne_id_callback(self, attr, old, new):
        print('switching from {} to {}'.format(old, new))
        self.ds_bgc.selected.indices = []
        self.ds_spec.selected.indices = []
        self.bgc_tsne_id = new
        # messy but it works and seems to be best way of doing this
        newdata = self.nh.bgc_data[new]
        self.bgc_datasource.data.update(
                                             x=newdata['x'],
                                             y=newdata['y'],
                                             name=newdata['name'],
                                             strain=newdata['strain'],
                                             gcf=newdata['gcf'],
                                             fill=newdata['fill'])

    def bokeh_layout(self):
        self.spec_div = Div(text="", sizing_mode='scale_height', name='spec_div')
        self.bgc_div = Div(text="", sizing_mode='scale_height', name='bgc_div')

        self.ds_bgc.selected.on_change('indices', self.bgc_selchanged)

        self.plot_toggles = CheckboxButtonGroup(active=[PLOT_CMAP, PLOT_SINGLETONS], labels=PLOT_TOGGLES, name='plot_toggles')
        self.plot_toggles.on_change('active', self.plot_toggles_callback)

        self.tsne_id_select = Select(title='BGC TSNE:', value=self.bgc_tsne_id, options=self.bgc_tsne_id_list, name='tsne_id_select')
        self.tsne_id_select.on_change('value', self.tsne_id_callback)

        self.scoring_mode_group = RadioGroup(labels=SCORING_MODES, active=SCO_MODE_BGC_SPEC, name='scoring_mode_group')
        self.scoring_mode_group.on_change('active', self.scoring_mode_callback)

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
        self.sliders = row(self.metcalf_percentile, self.likescore_cutoff, self.hg_prob, name='sliders')

        # for debug output etc 
        self.debug_div = Div(text="", name='debug_div')

        curdoc().add_root(self.plot_toggles)
        curdoc().add_root(self.tsne_id_select)
        curdoc().add_root(self.fig_spec)
        curdoc().add_root(self.fig_bgc)
        curdoc().add_root(self.spec_div)
        curdoc().add_root(self.bgc_div)
        curdoc().add_root(self.scoring_mode_group)
        curdoc().add_root(self.scoring_method_group)
        curdoc().add_root(self.sliders)
        curdoc().add_root(self.debug_div)

# server_lifecycle.py adds a .nh attr to the current Document instance, use that
# to access the already-created NPLinker instance plus the TSNE data
nb = NPLinkerBokeh(curdoc().nh)
nb.create_plots()
nb.bokeh_layout()

