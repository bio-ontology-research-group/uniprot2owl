from java.io import FileReader
import os
import time

from org.apache.jena.rdf.model import ModelFactory
from org.apache.jena.vocabulary import RDF
from org.semanticweb.elk.owlapi import ElkReasonerFactory
from org.semanticweb.owlapi.apibinding import OWLManager
from org.semanticweb.owlapi.model import IRI, AddImport, OWLClass, \
    OWLOntologyLoaderConfiguration, OWLOntologyManager
from org.semanticweb.owlapi.model import OWLObject
from org.semanticweb.owlapi.reasoner import ConsoleProgressMonitor, \
    SimpleConfiguration
from org.semanticweb.owlapi.vocab import OWLRDFVocabulary

from Queue import Queue
import collections
from threading import Thread
import traceback


# Number of threads to use
numthreads = 48

# Choose file directories here
input_directory = "/home/mencella/borg/swissprot_ttls/"
output_directory = "/home/mencella/borg/swissprot.owl"
gotaxon = "/home/mencella/borg/gotaxon.owl"

up = "http://purl.uniprot.org/core/"

manager = OWLManager.createOWLOntologyManager()
factory = manager.getOWLDataFactory()
print "Loading ontology..."
ontology = manager.createOntology(IRI.create("http://aber-owl.net/uniprot.owl"))
# ontology = manager.loadOntologyFromOntologyDocument(IRI.create("file:" + gotaxon))
print "... ontology loaded."
onturi = "http://aber-owl.net/uniprot.owl#"

# Imports
config = OWLOntologyLoaderConfiguration()
 
importDeclaraton = factory.getOWLImportsDeclaration(IRI.create("http://purl.obolibrary.org/obo/go.owl"));
manager.applyChange(AddImport(ontology, importDeclaraton))
manager.makeLoadImportRequest(importDeclaraton, config)
 
importDeclaraton = factory.getOWLImportsDeclaration(IRI.create("http://purl.obolibrary.org/obo/ncbitaxon.owl"));
manager.applyChange(AddImport(ontology, importDeclaraton))
manager.makeLoadImportRequest(importDeclaraton, config)
 
# manager.saveOntology(ontology, IRI.create("file:" + "/home/mencella/borg/gotaxon.owl"))

def create_relation(s):
    if s == "part-of":
        istring = "http://purl.obolibrary.org/obo/BFO_0000050"
    elif s == "has-part":
        istring = "http://purl.obolibrary.org/obo/BFO_0000051"
    else:
        istring = "http://phenomebrowser.net/#" + s
    return factory.getOWLObjectProperty(IRI.create(istring))
    
def create_class(s):
    return factory.getOWLClass(IRI.create(s))
    
def add_anno(resource, prop, cont):
    anno = factory.getOWLAnnotation(factory.getOWLAnnotationProperty(prop.getIRI()), factory.getOWLLiteral(cont))
    axiom = factory.getOWLAnnotationAssertionAxiom(resource.getIRI(), anno)
    manager.addAxiom(ontology, axiom)
    
# Subclasses of 'cellular location'
progressMonitor = ConsoleProgressMonitor();
config = SimpleConfiguration(progressMonitor);
reasoner = ElkReasonerFactory().createReasoner(ontology, config);
loc_cls = create_class("http://purl.obolibrary.org/obo/GO_0005575")
loc_nodeset = reasoner.getSubClasses(loc_cls, False).getFlattened()
print loc_nodeset
    
genericProteinNames = dict()
proteinCounter = 0

# Initiate a queue of all ttl files to be read.
queue = Queue()
for tfile in os.listdir(input_directory):
    queue.put(tfile)
print "Queue built. There are %d files to be read from %s." % (queue._qsize(), input_directory)

def readFiles(i, q):
    while True:
        tfile = q.get()
        
        size = q._qsize()
        if size % 1000 == 0:
            print "%d entries left in queue. " % size, time.strftime("%c")
            
        rdfModel = ModelFactory.createDefaultModel()
        try:
            rdfModel.read(FileReader(input_directory + tfile), "http://foobar#", "TURTLE")
        except:
            print "Error with file: ", tfile
            q.task_done()
            #traceback.print_exc()
            continue
        
        iri_iter = rdfModel.listStatements(None, RDF.type, rdfModel.createResource(up + "Protein"))
        
        while iri_iter.hasNext(): # iterate over all Protein iri's in file
            iris = iri_iter.nextStatement().getSubject()
            iri = iris.toString()
            label = rdfModel.listStatements(iris, rdfModel.createProperty(up + "mnemonic"), None).nextStatement().getObject().toString()
            
            # subclass
            cl = create_class(iri)
            add_anno(cl, OWLRDFVocabulary.RDFS_LABEL, label)
            cls = create_class(iri + "_all")
            add_anno(cls, OWLRDFVocabulary.RDFS_LABEL, "All %s in the universe" % label)
            manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(cl, factory.getOWLObjectSomeValuesFrom(create_relation("member-of"), cls)))
            
            # generic names
            genericNames = []
            first = True
            for stmt0 in rdfModel.listStatements(iris, rdfModel.createProperty(up + "submittedName"), None):
                name = stmt0.getObject()
                for stmt1 in rdfModel.listStatements(name, RDF.type, rdfModel.createResource(up + "Structured_Name")):
                    for stmt2 in rdfModel.listStatements(stmt1.getSubject(), rdfModel.createProperty(up + "fullName"), None):
                        if first:
                            first = False
                            firstName = stmt2.getObject().toString()
                        genericNames.append(stmt2.getObject().toString())
                        
                if firstName not in genericProteinNames.keys():
                    ncl = create_class(onturi+"GENERICPROTEIN_%s" % proteinCounter)
                    ncla = create_class(onturi+"GENERICPROTEIN_%s_all" % proteinCounter)
                    proteinCounter += 1
                    for name in genericNames:
                        universe = "All %s in the universe" % name
                        add_anno(ncl, OWLRDFVocabulary.RDFS_LABEL, name) # add label
                        add_anno(ncla, OWLRDFVocabulary.RDFS_LABEL, universe) # add label
                        genericProteinNames[name] = ncl
                        genericProteinNames[universe] = ncla
            
                ncl = genericProteinNames[firstName]
                ncla = genericProteinNames["All %s in the universe" % firstName]
                manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(cl, ncl))
                manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(cls, ncla))
            
            # isoforms
            isoforms = collections.OrderedDict()
            for stmt in rdfModel.listStatements(iris, rdfModel.createProperty(up + "sequence"), None):
                isoforms[stmt.getObject().toString()] = 1
            for iso in isoforms.keys():
                icl = create_class(iso)
                add_anno(icl, OWLRDFVocabulary.RDFS_LABEL, "Isoform of %s" % label) # add label
                manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(icl, factory.getOWLObjectSomeValuesFrom(create_relation("isoform-of"), cl)))
                manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(icl, cl))
                
            # gofunctions
            gofunctions = collections.OrderedDict()
            
            for stmt in rdfModel.listStatements(rdfModel.createResource(iri), rdfModel.createProperty(up+"classifiedWith"), None):
                if(stmt.getObject().toString().startswith("http://purl.obolibrary.org/obo/")):
                    gofunctions[stmt.getObject().toString()] = 1
            for fun in gofunctions.keys():
                # Decide whether this is a function/process or location
                if create_class(fun) in loc_nodeset: # location
                    relation = "located-in"
                else: # function/process
                    relation = "participates-in"
                manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(cls, factory.getOWLObjectSomeValuesFrom(create_relation("has-member"), 
                                                                        factory.getOWLObjectSomeValuesFrom(create_relation(relation), create_class(fun)))))
                
            organism = rdfModel.listStatements(rdfModel.createResource(iri), rdfModel.createProperty(up+"organism"), None).nextStatement().getObject().toString()
            prefix = "http://purl.obolibrary.org/obo/NCBITaxon_"
            slash = organism.rfind('/')
            organism = prefix + organism[slash+1:]
            manager.addAxiom(ontology, factory.getOWLSubClassOfAxiom(cl, factory.getOWLObjectSomeValuesFrom(create_relation("created-in-organism"), create_class(organism))))

        #signals to queue job is done
        q.task_done()


for i in range(numthreads):
    print "Thread %d initiated" % (i+1)
    t = Thread(target=readFiles, args=(i, queue))
    t.setDaemon(True)
    t.start()

start = time.clock()
print "Started at: ", time.strftime("%c")
queue.join()
print "Ended at: ", time.strftime("%c")
end = time.clock()

print "Conversion completed after %d seconds." % (end - start)
print "Number of axioms: ", ontology.getLogicalAxiomCount()

manager.saveOntology(ontology, IRI.create("file:" + output_directory))
print "Ontology saved as %s" % output_directory
