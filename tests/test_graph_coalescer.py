import pytest
import src.graph_coalescence.graph_coalescer as gc
import src.single_node_coalescer as snc
from src.components import Opportunity,Answer
import os,json


#Failing due to RK KG problems.  Once HGNC FAMILY is fixed, turn this back on.  CB May 6, 2020
def test_graph_coalescer():
    curies = [ 'NCBIGene:106632262', 'NCBIGene:106632263','NCBIGene:106632261' ]
    opportunity = Opportunity('hash',('qg_0','gene'),curies,[0,1,2])
    opportunities=[opportunity]
    patches = gc.coalesce_by_graph(opportunities)
    assert len(patches) == 1
    #patch = [qg_id that is being replaced, curies (kg_ids) in the new combined set, props for the new curies, answers being collapsed]
    p = patches[0]
    assert p.qg_id == 'qg_0'
    assert len(p.set_curies) == 3 # 3 of the 3 curies are subclasses of the output
    assert p.new_props['coalescence_method'] == 'graph_enrichment'
    assert p.new_props['p_value'] < 1e-10
    assert isinstance(p.new_props['enriched_nodes'],list)
    assert len(p.added_nodes)==1

def test_gene_taxon():
    """Do a bunch of genes return in taxon -> human?  That's a dumb result, but it is happening"""
    curies = [ "NCBIGene:7200", "NCBIGene:5805", "NCBIGene:7124", "NCBIGene:12", "NCBIGene:83985", "NCBIGene:54478",
               "NCBIGene:1401", "NCBIGene:160622", "NCBIGene:11113", "NCBIGene:790", "NCBIGene:10439", "NCBIGene:6651",
               "NCBIGene:640", "NCBIGene:4141", "NCBIGene:1677", "NCBIGene:5914", "NCBIGene:9104", "NCBIGene:5447",
               "NCBIGene:85480", "NCBIGene:1956", "NCBIGene:9457", "NCBIGene:10129", "NCBIGene:462", "NCBIGene:5978",
               "NCBIGene:5741", "NCBIGene:65057", "NCBIGene:920", "NCBIGene:50943", "NCBIGene:4810", "NCBIGene:196",
               "NCBIGene:828", "NCBIGene:1431", "NCBIGene:355", "NCBIGene:25759", "NCBIGene:4149", "NCBIGene:367",
               "NCBIGene:56005", "NCBIGene:7099", "NCBIGene:54440", "NCBIGene:1825", "NCBIGene:1674", "NCBIGene:6430",
               "NCBIGene:7054", "NCBIGene:8161", "NCBIGene:124872", "NCBIGene:847", "NCBIGene:80271", "NCBIGene:1636",
               "NCBIGene:3266", "NCBIGene:26762", "NCBIGene:1611", "NCBIGene:637", "NCBIGene:572", "NCBIGene:1154",
               "NCBIGene:2232", "NCBIGene:4056", "NCBIGene:9536", "NCBIGene:326", "NCBIGene:958", "NCBIGene:22858",
               "NCBIGene:6010", "NCBIGene:10603", "NCBIGene:2022", "NCBIGene:11092", "NCBIGene:6476", "NCBIGene:5742",
               "NCBIGene:126669", "NCBIGene:2056", "NCBIGene:3630", "NCBIGene:8685", "NCBIGene:5743", "NCBIGene:55806",
               "NCBIGene:4118", "NCBIGene:9235", "NCBIGene:831", "NCBIGene:929", "NCBIGene:9146", "NCBIGene:7450",
               "NCBIGene:5091", "NCBIGene:240", "NCBIGene:51324", "NCBIGene:4712" ]
    nc = len(curies)
    opportunity = Opportunity('hash',('qg_0','gene'),curies,range(nc))
    opportunities=[opportunity]
    patches = gc.coalesce_by_graph(opportunities)
    assert len(patches) == 1
    #patch = [qg_id that is being replaced, curies (kg_ids) in the new combined set, props for the new curies, answers being collapsed]
    p = patches[0]
    assert p.qg_id == 'qg_0'
    assert len(p.set_curies) == 3 # 3 of the 3 curies are subclasses of the output
    assert p.new_props['coalescence_method'] == 'graph_enrichment'
    assert p.new_props['p_value'] < 1e-10
    assert isinstance(p.new_props['enriched_nodes'],list)
    assert len(p.added_nodes)==1

def test_graph_coalescer_perf_test():
    from src.single_node_coalescer import coalesce
    import datetime

    # get a timestamp
    t1 = datetime.datetime.now()

    # get the path to the test file
    test_filename = os.path.join(os.path.abspath(os.path.dirname(__file__)),'EdgeIDAsStrAndPerfTest.json')

    # open the file and load it
    with open(test_filename,'r') as tf:
        incoming = json.load(tf)
        incoming = incoming['message']

    # call function that does property coalesce
    coalesced = coalesce(incoming, method='graph')

    # get the amount of time it took
    diff = datetime.datetime.now() - t1

    # it should be less than this
    assert(diff.seconds < 60)

    # loop through the query_graph return and insure that edge ids are strs
    for n in coalesced['query_graph']['edges']:
        assert(isinstance(n['id'], str))

    # loop through the knowledge_graph return and insure that edge ids are strs
    for n in coalesced['knowledge_graph']['edges']:
        assert(isinstance(n['id'], str))

    # loop through the knowledge_graph return and insure that nodes have a type list
    for n in coalesced['knowledge_graph']['nodes']:
        assert(isinstance(n['type'], list))

def test_graph_coalesce():
    """Make sure that results are well formed."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    testfilename = os.path.join(dir_path, 'famcov_new.json')
    with open(testfilename, 'r') as tf:
        answerset = json.load(tf)
        answerset = answerset['message']
    newset = snc.coalesce(answerset, method='graph')
    kgnodes = set([n['id'] for n in newset['knowledge_graph']['nodes']])
    for r in newset['results']:
        nbs = r['node_bindings']
        extra = False
        for nb in nbs:
            if nb['qg_id'].startswith('extra'):
                extra = True
            #Every node binding should be found somewhere in the kg nodes
            for kgid in nb['kg_id']:
                assert kgid in kgnodes
        assert extra

def test_graph_coalesce_strider():
    """Make sure that results are well formed."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    testfilename = os.path.join(dir_path, 'strider_relay_mouse.json')
    with open(testfilename, 'r') as tf:
        answerset = json.load(tf)
        answerset = answerset['message']
    newset = snc.coalesce(answerset, method='graph')
    for r in newset['results']:
        nbs = r['node_bindings']
        extra = False
        for nb in nbs:
            if nb['qg_id'].startswith('extra'):
                extra = True
        assert extra

def test_missing_node_norm():
    from src.single_node_coalescer import coalesce
    import datetime

    # get a timestamp
    t1 = datetime.datetime.now()

    # get the path to the test file
    test_filename = os.path.join(os.path.abspath(os.path.dirname(__file__)),'graph_named_thing_issue.json')

    # open the file and load it
    with open(test_filename,'r') as tf:
        incoming = json.load(tf)
        incoming = incoming['message']

    # call function that does property coalesce
    coalesced = coalesce(incoming, method='graph')

    # get the amount of time it took
    diff = datetime.datetime.now() - t1

    # it should be less than this
    assert(diff.seconds < 60)

    # loop through the query_graph return and insure that edge ids are strs
    for n in coalesced['query_graph']['edges']:
        assert(isinstance(n['id'], str))

    # loop through the knowledge_graph return and insure that edge ids are strs
    for n in coalesced['knowledge_graph']['edges']:
        assert(isinstance(n['id'], str))

def test_gouper():
    x = 'abcdefghi'
    n = 0
    for group in gc.grouper(3,x):
        x = group
        n += 1
    assert n == 3
    assert x == ('g','h','i')

def test_gouper_keys():
    d = {x:x for x in 'abcdefg'}
    n = 0
    for group in gc.grouper(3, d.keys()):
        x = group
        n += 1
    assert n == 3
    assert x == ('g',)
