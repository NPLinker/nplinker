from bokeh.layouts import row
from bokeh.models import Div

# each of these is a 2-tuple, (title, text)
TOOLTIP_TABLES_FILTERING = (
    'Initial Metcalf filtering',
    """NPLinker runs Metcalf scoring on the dataset when initially loaded, and removes all objects scoring under the indicated threshold before loading the remainder into the tables. This is done to reduce the number of uninteresting objects to search through, and to make the application more responsive. You can adjust this value in the configuration file if required, see the user guide for details"""
)

TOOLTIP_SPECTRA_SCORING = (
    'Run scoring on spectra',
    """Clicking this button runs the currently enabled scoring method(s) on the set of visible spectra in the table below. Before using this, you need to make at least one selection from the Molecular Family or Spectrum tables. Results are displayed at the bottom of the page"""
)

TOOLTIP_BGC_SCORING = (
    'Run scoring on BGCs',
    """Clicking this button runs the currently enabled scoring method(s) on the set of visible BGCs in the table below. Before using this, you need to make at least one selection from the GCF or BGC tables. Results are displayed at the bottom of the page"""
)

TOOLTIP_CLEAR_SELECTIONS = (
    'Reset tables',
    """Click this button to clear all selections from all tables and go back to the original state"""
)

TOOLTIP_DOWNLOAD_CSV_MOLFAM = (
    'Download CSV data',
    """Download a CSV file containing the information currently visible in the MolFam table"""
)

TOOLTIP_DOWNLOAD_CSV_SPECTRA = (
    'Download CSV data',
    """Download a CSV file containing the information currently visible in the Spectra table"""
)

TOOLTIP_DOWNLOAD_CSV_BGC = (
    'Download CSV data',
    """Download a CSV file containing the information currently visible in the BGC table"""
)

TOOLTIP_DOWNLOAD_CSV_GCF = (
    'Download CSV data',
    """Download a CSV file containing the information currently visible in the GCF table"""
)

TOOLTIP_RESULTS_AREA = (
    'Scoring results',
    """Results from clicking either of the two scoring buttons above the tables will be displayed here"""
)

TOOLTIP_HYBRID_SCORING = (
    'Hybrid BGC handling',
    """When you run scoring with a selection of GCFs that include one or more hybrid BGCs with this option disabled, NPLinker will look for links among ALL parent GCFs of each hybrid BGC. This may mean you find GCFs in the results that you didn't originally select. For example, if GCF A contains a hybrid BGC with parents GCF A and GCF B, you may find GCF B in the results as it will be included in the search for links through the hybrid BGC. If you do NOT want this behaviour and instead want to limit NPLinker to searching for links only in the GCF(s) selected, enable this option. This will force NPLinker to ignore other parent GCFs associated with each hybrid BGC. If your input is a GCF with no hybrid BGCs, this option has no effect."""
)

_popover_text = """
    <a tabindex="0" class="btn btn-sm btn-secondary popover-dismiss" role="button" data-toggle="popover" data-trigger="focus" title="{}" data-content="{}">?</a>
"""
def create_popover(title, text):
    return _popover_text.format(title, text)

def wrap_popover(host, pop_div, name, sizing_mode='scale_width'):
    return row(children=[host, Div(text=pop_div)], name=name, sizing_mode=sizing_mode)
