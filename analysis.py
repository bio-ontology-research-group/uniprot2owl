import time

from org.semanticweb.elk.owlapi import ElkReasonerFactory
from org.semanticweb.owlapi.apibinding import OWLManager
from org.semanticweb.owlapi.model import IRI, AxiomType
from org.semanticweb.owlapi.reasoner import SimpleConfiguration, \
    ConsoleProgressMonitor, InferenceType

import random
import traceback


PREFIX_MAP = {
    "http://www.geneontology.org/formats/oboInOwl#" : "oboInOwl:",
    "http://purl.org/dc/elements/1.1/" : "dc:",
    "http://protege.stanford.edu/plugins/owl/protege#" : "protege:",
    "http://purl.org/dc/terms/" : "dc:",
    "http://purl.org/dc/elements/1.1/" : "dc:",
    "http://purl.org/dc/terms/" : "dc:",
    "http://purl.org/vocab/vann/" : "vann:",
    "http://purl.org/spar/cito/" : "cito:",
    "http://purl.obolibrary.org/obo/" : "obo:",
    "http://www.w3.org/2004/02/skos/core" : "skos:",
    "http://semanticscience.org/resource/" : "sio:"
}

def prefixUrls(s):
    for prefix in PREFIX_MAP.keys():
        if(s.startsWith(prefix)):
            s = s.replace(prefix, PREFIX_MAP[prefix])
    return s


file = "/home/mencella/borg/owl_files/swissprot.owl"
gofile = "/home/mencella/borg/owl_files/go.owl"
ncbifile = "/home/mencella/borg/owl_files/ncbitaxon.owl"


manager = OWLManager.createOWLOntologyManager()
factory = manager.getOWLDataFactory()
print "Loading OWL files..."
ontology = manager.loadOntologyFromOntologyDocument(IRI.create("file:" + file))
go = manager.loadOntologyFromOntologyDocument(IRI.create("file:" + gofile))
ncbi = manager.loadOntologyFromOntologyDocument(IRI.create("file:" + ncbifile))
print "... OWL files loaded."


def create_relation(s):
    return factory.getOWLObjectProperty(IRI.create("http://phenomebrowser.net/#" + s))

def create_class(s):
    return factory.getOWLClass(IRI.create(s))


stats = {
        'rBoxAxiomCount': 0,
        'tBoxAxiomCount': 0,
        'totalAxiomCount': 0,
        'unsatisfiableClassesCount':0,
        'logicalAxiomCount': 0,
        'complexity': 0,
        'annotations': dict(),
        'classCount': ontology.getClassesInSignature(True).size(),
        'loaded': (ontology in manager.ontologies),
#         'consistent': True #(ontology in manager.queryEngines)
        }

# unsatis = manager.runQuery("<http://www.w3.org/2002/07/owl#Nothing>", "equivalent", ontology, -1, False, False)
# stats.unsatisfiableClassesCount = unsatis.size()

for axtype in AxiomType.TBoxAxiomTypes:
    for ax in ontology.getAxioms(axtype, True):
        stats['totalAxiomCount'] += 1
        stats['tBoxAxiomCount'] += 1

for axtype in AxiomType.RBoxAxiomTypes:
    for ax in ontology.getAxioms(axtype, True):
        stats['totalAxiomCount'] += 1 
        stats['rBoxAxiomCount'] += 1 
  
stats['complexity'] = stats['totalAxiomCount'] / stats['classCount']
stats['logicalAxiomCount'] = ontology.getLogicalAxiomCount()

for anno in ontology.getAnnotations():
    try:
        prop = anno.getProperty().toString().replaceAll("<","").replaceAll(">","")
        prop = prefixUrls(prop)
        val = anno.getValue().getLiteral().toString()
        stats['annotations'][prop] = val
    except:
        traceback.print_exc()

print stats

progressMonitor = ConsoleProgressMonitor()
config = SimpleConfiguration(progressMonitor)
reasoner = ElkReasonerFactory().createReasoner(ontology)
thing = factory.getOWLThing()

# Measure classification time
total_time = 0
repeats = 10
for i in range(repeats):
    start = time.clock()
    print "Start timer..."
    reasoner = ElkReasonerFactory().createReasoner(ontology)
    reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)
    end = time.clock()
    total_time += (end - start)
    print "Classification time: %d seconds." % (end - start)
print "Average classification time: %d seconds" % (total_time / repeats)


#### Query 1

total_time1, total_time2 = 0, 0
repeats = 1000

go_reasoner = ElkReasonerFactory().createReasoner(go)
go_top_set = go_reasoner.getSubClasses(thing, False).getFlattened()
go_top_list = []
for x in go_top_set:
    go_top_list.append(x)
print "%d GO classes" % len(go_top_list)

ncbi_reasoner = ElkReasonerFactory().createReasoner(ncbi)
ncbi_top_set = ncbi_reasoner.getSubClasses(thing, False).getFlattened()
ncbi_top_list = []
for x in ncbi_top_set:
    ncbi_top_list.append(x)
print "%d NCBI classes" % len(ncbi_top_list)

# Cellular location
loc_cls = create_class("http://purl.obolibrary.org/obo/GO_0005575")
loc_nodeset = reasoner.getSubClasses(loc_cls, False).getFlattened()

print "Testing queries..."
for i in range(repeats):
    print "%d queries left" % (repeats - i)
    go_cl = random.choice(go_top_list)
    if go_cl in loc_nodeset:   # function/process
        relation = "located-in"
    else:
        relation = "participates-in"
    
    temp1 = factory.getOWLObjectSomeValuesFrom(create_relation(relation), go_cl)
    temp2 = factory.getOWLObjectSomeValuesFrom(create_relation("has-member"), temp1)
    query_class = factory.getOWLObjectSomeValuesFrom(create_relation("member-of"), temp2)
    start = time.clock()
    reasoner.getSubClasses(query_class, False)
    end = time.clock()
    total_time1 += (end - start)
    
    ncbi_cl = random.choice(ncbi_top_list)
    temp = factory.getOWLObjectSomeValuesFrom(create_relation("created-in"), ncbi_cl)
    query2 = factory.getOWLObjectIntersectionOf(query_class, temp)
    start = time.clock()
    reasoner.getSubClasses(query2, False)
    end = time.clock()
    total_time2 += (end - start)
    
print "Average query 1 time: %f seconds" % ((total_time1 + 0.0) / repeats)
print "Average query 2 time: %f seconds" % ((total_time2 + 0.0) / repeats)


print "Program terminated."
