import pytest
import json
import os
import src.single_node_coalescer as snc
import src.ontology_coalescence.ontology_coalescer as oc
from src.components import Opportunity


def test_big_ontology():
    fn = 'bigger.json'
    testfilename = os.path.join(os.path.abspath(os.path.dirname(__file__)),fn)
    with open(testfilename,'r') as tf:
        answerset = json.load(tf)
    newset = snc.coalesce(answerset,method='ontology')
    rs = newset['results']
    print(len(rs))

def test_big_graphbased():
    fn = 'bigger.json'
    testfilename = os.path.join(os.path.abspath(os.path.dirname(__file__)),fn)
    with open(testfilename,'r') as tf:
        answerset = json.load(tf)
    newset = snc.coalesce(answerset,method='graph')
    rs = newset['results']
    print(len(rs))
