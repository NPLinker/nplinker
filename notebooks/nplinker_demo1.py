#!/usr/bin/env python

# In[1]:


# if running from clone of the git repo
# sys.path.append('../src')

# import the main NPLinker class. normally this all that's required to work
# with NPLinker in a notebook environment
from nplinker.nplinker import NPLinker


# In[2]:


# the standard method of loading a dataset configuration is to pass the filename
# of a TOML configuration file to the NPLinker constructor.
npl = NPLinker('./nplinker_demo1.toml')
# loading the actual data files can take some time depending on the dataset,
# so this is done separately by calling the load_data method.
#
# During the loading process, logging messages will be printed to stdout. This
# can be useful for debugging problems with files not being discovered or parsed
# correctly. You can control the verbosity of these messages in the configuration
# file if required, and/or redirect them to a file instead of stdout.
npl.load_data()


# In[ ]:


# Basic functionality
# ===================
#
# Once you have an NPLinker object with all data loaded, there are a collection of simple
# methods and properties you can use to access objects and metadata. Some examples are
# given below, see https://nplinker.readthedocs.io/en/latest/ for a complete API description.

# configuration/dataset metadata
# - a copy of the configuration as parsed from the .toml file (dict)
print(npl.config)
# - the path to the directory where various nplinker data files are located (e.g. the
#   default configuration file template) (str)
print(npl.data_dir)
# - a dataset ID, derived from the path for local datasets or the paired platform ID
#   for datasets loaded from that source (str)
print(npl.dataset_id)
# - the root directory for the current dataset (str)
print(npl.root_dir)

# objects
# - you can directly access lists of each of the 4 object types:
print('BGCs:', len(npl.bgcs))
print('GCFs:', len(npl.gcfs)) # contains GCF objects
print('Spectra:', len(npl.spectra)) # contains Spectrum objects
print('Molecular Families:', len(npl.molfams)) # contains MolecularFamily objects


# In[ ]:


# Scoring functionality - part 1
# ==============================
# (again see https://nplinker.readthedocs.io/en/latest/ for API documentation)

# NPLinker provides a set of scoring methods that can be used individually or
# in combination to find interesting links in the current dataset. To get a
# get a list of the names of the available scoring methods:
print('Available scoring methods:')
for m in npl.scoring_methods:
    print(f' - {m}')

# The first step in running a scoring operation is to get an instance of the
# method(s) you want to use by calling scoring_method():
mc = npl.scoring_method('metcalf')

# Now mc is an instance of the class that implements Metcalf scoring. Once
# you have such an instance, you may change any of the parameters it exposes.
# In the case of Metcalf scoring, the following parameters are currently exposed:
# - cutoff (float): the scoring threshold. Links with scores less than this are excluded
# - standardised (bool): set to True to use standardised scores (default), False for regular
mc.cutoff = 3.5
mc.standardised = True


# In[ ]:


# Scoring functionality - part 2
# ==============================

# After creating and optionally configuring a scoring method, you need to call
# get_links() to perform the operation on a selected set of objects. This method
# takes 2-3 parameters, the third being optional:
#  - a list of objects to find links from (or a list of lists of objects)
#  - a list of scoring methods, or a single method as shorthand for a 1-element list
#  - (optional) a boolean indicating if results from multiple methods should be
#     ANDed together to produce the final results. If set to False, the results will
#     contain links found by any method rather than all methods.
#
# This first example shows the simplest case: 1 set of objects and 1 scoring method.
# If the and_mode parameter is not given it defaults to True, but the value doesn't
# matter here because only one method is being used.
results = npl.get_links(npl.gcfs[:10], mc, and_mode=True)

# get_links returns an instance of a class called LinkCollection. This provides a wrapper
# around the results of the scoring operation and has various useful properties/methods:
#
# - len(results) or .source_count will tell you how many of the input_objects were found to have links
print(f'Number of results: {len(results)}')
# - .sources is a list of those objects
objects_with_links = results.sources
# - .links is a dict with structure {input_object: {linked_object: ObjectLink}}
objects_and_link_info = results.links
# - .get_all_targets() will return a flat list of *all* the linked objects (for all sources)
all_targets = results.get_all_targets()
# - .methods is a list of the scoring methods passed to get_links
methods = results.methods


# In[ ]:


# Scoring functionality - part 3
# ==============================
#
# The link data inside the LinkCollection object is itself stored in ObjectLink objects.
# Each instance of an ObjectLink represents a link between a given pair of objects as
# determined by 1 or more scoring methods.
#
# ObjectLinks have the following basic attributes:
# - .source: the input object provided to the method
# - .target: the linked object
# - .methods: a list of the methods that found this link
# - .shared_strains: a list of Strain objects (possibly empty) shared between .source and .target
# - .data(<method_object>): return the output of <method_object> for this link (e.g. any score values)
#
# You can also retrieve any method-specific info for a link by subscripting these objects with
# the appropriate method object, e.g. metcalf_link_data = object_link[mc]

# This shows how to iterate over the link information from result.links. In the body of the loop
# <obj> will be one of  the original objects supplied to get_links and <result> will be a dict
# with structure {linked_object: ObjectLink} (indicating <obj> is linked to <linked_object> according to
# the information stored in the ObjectLink)
for obj, result in results.links.items():
    # display the object, the number of links it has, and the number of methods that were used to get them
    print(f'Results for object: {obj}, {len(result)} total links, {results.method_count} methods used')

    # sorting is method-dependent since they might have very different "scores", so you should
    # use the original object to do this. For Metcalf scoring, this will return the ObjectLinks sorted
    # by their Metcalf scores.
    sorted_links = results.get_sorted_links(mc, obj)
    # or if you wanted them in the reverse order:
    # sorted_links = results.get_sorted_links(mc, obj, reverse=True)

    # Now display some link information for each link associated with <obj>.
    # link_data[<method_object>] will return the per-link data generated by that
    # method. Here the metcalf method simply returns the link score as a floating point value,
    # but other methods may return more complex objects.
    #
    # Each scoring method also has a format_data method which should provide a relatively short
    # human-readable summary of the data, as a quick way to print and examine results.
    for link_data in sorted_links:
        print('  --> [{}] {} | {} | shared strains = {}'.format(','.join(method.name for method in link_data.methods),
                                                                link_data.target,
                                                                mc.format_data(link_data[mc]),
                                                                len(link_data.shared_strains)))

    # alternatively, if you don't care about ordering, you can just iterate directly over the
    # linked objects like this:
    # for link_target, link_data in result.items():
    #    print(link_target, link_data)



# In[ ]:


# Scoring functionality - part 4
# ==============================
#
# The LinkCollection object supports performing various types of filtering on the original set of
# results contained within it:
# - .filter_no_shared_strains(): remove any links where the linked objects do not share strains
# - .filter_sources(callable), .filter_targets(callable), .filter_links(callable): each of these
#     simply execute callable(object) and filter out objects for which the return value is False/0.
#     The <objects> in each case are respectively: the original input objects (sources),
#     their linked objects (targets), and the ObjectLink objects (links).
#
# NOTE:
# - these methods all modify the original LinkCollection in-place
# - they will automatically remove any original results for which no links exist after filtering. For
#    example, if there is a source object which starts off with 2 links, but has 0 after a filter is
#    run, this object will not appear in the LinkCollection afterwards.
#
# Examples:
# - exclude any sources for which an arbitrary function is false (sources are GCFs in this example)
results.filter_sources(lambda gcf: gcf.id % 2 == 0)
# - exclude any linked objects for which an arbitrary function is false (targets are Spectrum objects here)
results.filter_targets(lambda spec: spec.id % 1 == 0)
# - exclude any links for which an arbitrary function is false (<link> is an ObjectLink)
results.filter_links(lambda link: link[mc] > 3.6)


# In[ ]:


# Scoring functionality - part 5
# ==============================
#
# The get_links method can be passed more complex parameters types than the above example which
# used a flat list of input objects and a single scoring method instance.

ts = npl.scoring_method('testscore') # copy of Metcalf method, only for debug use

# You can use the same set of objects with two different methods, and AND the results
# together so that objects will only be returned which have links according to
# BOTH of the supplied methods (if you provide 2 or more scoring methods but only a single
# set of objects, that set will be used as input to every method).
results = npl.get_links(npl.gcfs[:10], [mc, ts], and_mode=True)
