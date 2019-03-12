import os, sys
from random import random

from bokeh.models.widgets import RadioGroup, Slider, Div, PreText
from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import Button
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.events import LODEnd, LODStart

from nplinker import Spectrum

# TOOLS="crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
TOOLS="crosshair,pan,wheel_zoom,box_zoom,reset,tap,save,box_select"

PW, PH = 500, 500

def xrange_cb(attr, old, new):
    print('X', attr, old, new)

def yrange_cb(attr, old, new):
    print('Y', attr, old, new)

class Range(object):

    def __init__(self, s=None, e=None):
        self.s = s
        self.e = e

    def update(self, s, e):
        self.s = s
        self.e = e

    def range(self):
        if self.s is None or self.e is None:
            return 0

        return self.e - self.s

    def __repr__(self):
        return 'Range(s={}, e={}, r={})'.format(self.s, self.e, self.range())

class ZoomHandler(object):
    """
    This class deals with scaling the circle glyphs on the scatterplots so that
    they stay roughly the same size during zooming. 
    """

    def __init__(self, plot, renderer, threshold=0.05):
        self.x = Range()
        self.ox = Range()
        self.y = Range()
        self.oy = Range()
        self.threshold = threshold
        self.lod_mode = False
        self.renderer = renderer
        self.ratio = 1

        # want to listen for changes to the start/end of the x/y ranges displayed by the plot...
        plot.x_range.on_change('start', self.x_cb)
        plot.x_range.on_change('end', self.x_cb)
        plot.y_range.on_change('start', self.y_cb)
        plot.y_range.on_change('end', self.y_cb)

        # ... and also entering/leaving "level of detail" mode, a lower-detail mode activated
        # during panning/zooming etc. this allows changes producing by panning which we don't
        # care about because it doesn't alter the zoom level
        plot.on_event(LODStart, self.lod_event)
        plot.on_event(LODEnd, self.lod_event)
        
        self.plot = plot

    def xchanged(self):
        """
        Return True if the x range has changed
        """
        xd = abs(self.x.range() - self.ox.range())
        if xd > (self.threshold * self.ox.range()):
            return True

        return False

    def ychanged(self):
        """
        Return True if the y range has changed
        """
        yd = abs(self.y.range() - self.oy.range())
        if yd > (self.threshold * self.oy.range()):
            return True

        return False

    def x_cb(self, attr, old, new):
        """ 
        Callback for changes to the x range. attr = 'start'/'end', 'old' and 'new'
        give the previous and current values.
        """
        if attr == 'start':
            self.ox.s = old
            self.x.s = new
        else:
            self.ox.e = old
            self.x.e = new

        if old is None and attr == 'end':
            self.ratio = 1/150.0 
            self.renderer.glyph.radius = self.ratio * abs(self.plot.x_range.end - self.plot.x_range.start)
            print('Radius is now: {}'.format(self.renderer.glyph.radius))

        if attr == 'end' and not self.lod_mode and self.xchanged():
            print('CHANGE x range: {} v {}, {}'.format(self.x.range(), self.ox.range(), abs(self.x.range() - self.ox.range())))
            self.renderer.glyph.radius = self.ratio * self.x.range()

    def y_cb(self, attr, old, new):
        """ 
        Callback for changes to the y range. attr = 'start'/'end', 'old' and 'new'
        give the previous and current values.
        """
        if attr == 'start':
            self.oy.s = old
            self.y.s = new
        else:
            self.oy.e = old
            self.y.e = new

        if attr == 'end' and not self.lod_mode and self.ychanged():
            print('CHANGE y range: {} v {}, {}'.format(self.y.range(), self.oy.range(), abs(self.y.range() - self.oy.range())))

    def lod_event(self, event):
        self.lod_mode = isinstance(event, LODStart)

class NPLinkerBokeh(object):

    def __init__(self, helper):
        self.nh = helper

        self.bgc_datasource = ColumnDataSource(data=self.nh.bgc_data)
        self.spec_datasource = ColumnDataSource(data=self.nh.spec_data)

        self.bgc_zoom = None
        self.spec_zoom = None

    @property
    def nplinker(self):
        return self.nh.nplinker

    @property
    def bgc_data(self):
        return self.nh.bgc_data

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

        # TODO select initial radius somehow

        # create the BGC figure and populate it from the datasource (note the dict keys are used in the 
        # call to circle to identify attributes)

        radius = 0.1

        f_bgc = figure(tools=TOOLS, title="BGCs (n={})".format(len(self.nh.bgc_data['x'])), sizing_mode='scale_width', name="fig_bgc")
        r_bgc = f_bgc.circle('x', 'y', source=self.bgc_datasource, 
                radius=radius,
                radius_dimension='max',
                fill_alpha=0.6, 
                line_color=None,
                selection_color='#ff0000')

        self.bgc_zoom = ZoomHandler(f_bgc, r_bgc)

        # customize hover tooltip to show BGC name and index
        # (there are some special builtin properties prefixed with $, and then 
        # properties prefixed with @ come from the data source)
        hover_bgc = HoverTool(tooltips=[('index', '$index'), ('BGC name', '@name'), ('GCF name', '@gcf')])
        f_bgc.add_tools(hover_bgc)

        # create the MolFam figure in the same way
        # f_spec = figure(tools=TOOLS, title="Spectra (n={})".format(len(self.nh.spec_data['x'])), plot_width=PW, plot_height=PH, name="fig_spec")
        f_spec = figure(tools=TOOLS, title="Spectra (n={})".format(len(self.nh.spec_data['x'])), sizing_mode='scale_width', name="fig_spec")
        r_spec = f_spec.circle('x', 'y', source=self.spec_datasource, 
                            fill_alpha=0.6, 
                            radius=radius,
                            radius_dimension='max',
                            fill_color='#449944', 
                            line_color=None,
                            selection_color='#ff7f00', 
                            selection_fill_alpha=0.9,
                            nonselection_fill_alpha=0.05, 
                            nonselection_fill_color='#333333', 
                            nonselection_line_color=None, 
                            nonselection_line_alpha=0)

        self.spec_zoom = ZoomHandler(f_spec, r_spec)

        hover_spec = HoverTool(tooltips=[('index', '$index'), ('Spectra name', '@name')])
        f_spec.add_tools(hover_spec)

        return f_bgc, r_bgc, f_spec, r_spec

# server_lifecycle.py adds a .nh attr to the current Document instance, use that
# to access the already-created NPLinker instance plus the TSNE data
nb = NPLinkerBokeh(curdoc().nh)
nplinker = nb.nh.nplinker

fig_bgc, ren_bgc, fig_spec, ren_spec = nb.create_plots()
ds_spec = ren_spec.data_source
ds_bgc = ren_bgc.data_source

def update_debug_output():
    debug_div.text = '<pre>{}</pre>'.format(str(nplinker.scoring))

# format params:
# - id for header, e.g spec_heading_1
# - colour for header
# - id for body, e.g. spec_body_1
# - title text
# - id for body, same as #2
# - parent id 
# - main text
tmpl = '''
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

def update_bgc_output():
    test = '<h4>{} selected BGCs</h4>'.format(len(ds_bgc.selected.indices))
    for i in ds_bgc.selected.indices[:1]: # TODO multi-obj
        print('bgc index {}'.format(i))
        bgc = results.bgc
        gcf = bgc.parent

        hdr_id = 'bgc_header_{}'.format(i)
        body_id = 'bgc_body_{}'.format(i)
        title = 'BGC #{}, name={}'.format(i, nb.bgc_data['name'][i])
        body = '<dl>'
        for attr in ['name', 'strain']:
            body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(results.bgc, attr))
        body += '<dt>Parent GCF:</dt><dd>{}</dd>'.format(gcf)
        body += '</dl>'

        test += tmpl.format(hdr_id, 'b4c4e8', body_id, title, body_id, 'accordionBGC', body)

    bgc_div.text = test

def update_spec_output():
    test = '<h4>{} selected spectra</h4>'.format(len(ds_spec.selected.indices))
    # for i in ds_spec.selected.indices:
    if results.links is not None:
        for i, spec_score in enumerate(results.links):
            # spec = nplinker.spectra[i]
            spec, score = spec_score
            print(spec, score)
            hdr_id = 'spec_header_{}'.format(i)
            body_id = 'spec_body_{}'.format(i)
            title = 'id={}, score={}'.format(spec.id, score)
            body = '<dl>'
            for attr in ['id', 'spectrum_id', 'family']:
                body += '<dt>{}:</dt><dd>{}</dd>'.format(attr, getattr(spec, attr))
            body += '<dt>shared strains</dt><dd>{} (with {})</dd>'.format(results.shared_strains[spec], results.gcf)
            body += '</dl>'
            test += tmpl.format(hdr_id, 'adeaad', body_id, title, body_id, 'accordionSpec', body)
    spec_div.text = test

class CurrentResults(object):
    
    def __init__(self):
        self.clear()

    def clear(self):
        self.bgc = None
        self.gcf = None
        self.links = []
        self.shared_strains = {}

    def update(self, bgc, gcf, links, shared_strains):
        self.bgc = bgc
        self.gcf = gcf
        self.links = links
        self.shared_strains = {}
        if shared_strains is not None:
            for x in links:
                self.shared_strains[x[0]] = shared_strains[(x[0], gcf)]

results = CurrentResults()

def get_links(bgc_index):
    print('Selected BGC is {}'.format(bgc_index))
    bgc = nplinker.lookup_bgc(nb.bgc_data['name'][bgc_index])
    gcf = bgc.parent
    print('Contains {} BGCs'.format(len(gcf.bgc_list)))
    # 3. use nplinker to get a list of links from this GCF to Spectra
    nplinker.clear_links() # TODO hack
    current_method = nplinker.scoring.enabled()[scoring_method_group.active]
    print(current_method)
    links = nplinker.get_links(gcf, scoring_method=current_method)
    if len(links) > 0:
        spectra_and_scores = nplinker.links_for_obj(gcf, current_method, type_=Spectrum)
        print('NUM SPECTRA == {}'.format(len(spectra_and_scores)))
        spectra_point_ids = [nb.spec_indices[spec[0].spectrum_id] for spec in spectra_and_scores]
        # 4. take the indexes of the corresponding spectra from the datasource and highlight on the other plot
        ds_spec.selected.indices = spectra_point_ids

        shared_strains = nplinker.get_common_strains([gcf], [spec[0] for spec in spectra_and_scores])
        results.update(bgc, gcf, spectra_and_scores, shared_strains)
    else:
        print('No links found!')
        ds_spec.selected.indices = []
        results.clear()
        results.update(bgc, gcf, None, None)
    
    update_bgc_output()
    update_spec_output()
    update_debug_output()

def bgc_selchanged(attr, old, new):
    """
    Callback for changes to the BGC datasource selected indices
    """
    print('selection changed, old={}, new={}'.format(old, new))

    # if the new selection contains any valid indices (it may be empty)
    if len(new) > 0:
        # pick the first one (TODO: support multiple selections)
        bgc_index = new[0] 
        # get links for the selected BGC and update the plots
        get_links(bgc_index)
    else:
        # if selection is now empty, clear the selection on the spectra plot too
        ds_spec.selected.indices = []


# set up a callback for changes to selected indices in the BGC datasource
ds_bgc.selected.on_change('indices', bgc_selchanged)
# TODO need same for spectra

def metcalf_percentile_callback(attr, old, new):
    nplinker.scoring.metcalf.sig_percentile = new
    if len(ds_bgc.selected.indices) > 0:
        get_links(ds_bgc.selected.indices[0])

    update_debug_output()

def likescore_cutoff_callback(attr, old, new):
    nplinker.scoring.likescore.cutoff = new / 100.0
    if len(ds_bgc.selected.indices) > 0:
        get_links(ds_bgc.selected.indices[0])

    update_debug_output()

def hg_prob_callback(attr, old, new):
    nplinker.scoring.hg.prob = new / 100.0
    if len(ds_bgc.selected.indices) > 0:
        get_links(ds_bgc.selected.indices[0])

    update_debug_output()

def scale_bgc_callback(attr, old, new):
    s = new / 100.0
    print('scale_bgc {}'.format(s))
    ds_bgc.data['radius'] = [s for i in range(len(ds_bgc.data['radius']))]

def scale_spec_callback(attr, old, new):
    s = new / 100.0
    print('scale_spec {}'.format(s))
    ds_spec.data['radius'] = [s for i in range(len(ds_spec.data['radius']))]

def scoring_method_callback(attr, old, new):
    if len(ds_bgc.selected.indices) > 0:
        get_links(ds_bgc.selected.indices[0])
    update_debug_output()

# TODO tomorrow:
# - displaying proper lists of objects
#   - use Div() objects, added as roots and embedded in the bootstrap template

# TODO fix these
# display selected object info
spec_div = Div(text="", sizing_mode='scale_height', name='spec_div')
bgc_div = Div(text="", sizing_mode='scale_height', name='bgc_div')

scale_spec = Slider(start=1, end=100, value=int(100 * ds_spec.data['radius'][0]), step=1, title='Spectra scatterplot scaling', name='scale_spec')
scale_bgc = Slider(start=1, end=100, value=int(100 * ds_bgc.data['radius'][0]), step=1, title='BGC scatterplot scaling', name='scale_bgc')
scale_spec.on_change('value', scale_spec_callback)
scale_bgc.on_change('value', scale_bgc_callback)

scoring_method_group = RadioGroup(labels=[m for m in nplinker.scoring.enabled_names()], active=0, name='scoring_method_group')
scoring_method_group.on_change('active', scoring_method_callback)
# scoring_select = RadioGroup(labels=[m for m in nplinker.scoring.all_names()], active=0)
# metcalf stuff
metcalf_percentile = Slider(start=70, end=100, value=nplinker.scoring.metcalf.sig_percentile, step=1, title='[metcalf] sig_percentile')
metcalf_percentile.on_change('value', metcalf_percentile_callback)
# likescore
likescore_cutoff = Slider(start=0, end=100, value=int(100 * nplinker.scoring.likescore.cutoff), step=1, title='[likescore] cutoff x 100 = ')
likescore_cutoff.on_change('value', likescore_cutoff_callback)
# hg
hg_prob = Slider(start=0, end=100, value=int(100 * nplinker.scoring.hg.prob), step=1, title='[hg] prob x 100 = ')
hg_prob.on_change('value', hg_prob_callback)
sliders = row(metcalf_percentile, likescore_cutoff, hg_prob, name='sliders')

# for debug output etc 
debug_div = Div(text="", name='debug_div')


curdoc().add_root(fig_spec)
curdoc().add_root(fig_bgc)
curdoc().add_root(spec_div)
curdoc().add_root(bgc_div)
curdoc().add_root(scale_spec)
curdoc().add_root(scale_bgc)
curdoc().add_root(scoring_method_group)
curdoc().add_root(sliders)
curdoc().add_root(debug_div)
