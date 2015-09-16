@Grab(group='org.apache.jena', module='jena-core', version='3.0.0')
@Grab(group='org.semanticweb.elk', module='elk-owlapi', version='0.4.2')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.0.2')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.0.2')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.0.2')

import org.apache.jena.vocabulary.*
import org.apache.jena.rdfxml.xmlinput.*
import org.apache.jena.rdf.model.*
import org.semanticweb.owlapi.io.* 
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary

def up = "http://purl.uniprot.org/core/"
def isoform="http://purl.uniprot.org/isoforms"
def go = "http://purl.obolibrary.org/obo/"
def base = "http://purl.uniprot.org/uniprot"
def taxon = "http://purl.uniprot.org/taxonomy/"

RDFReader mr = new JenaReader()

OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
OWLDataFactory fac = manager.getOWLDataFactory()
def factory = fac

OWLOntology ont = manager.createOntology(IRI.create("http://aber-owl.net/uniprot.owl"))
def onturi = "http://aber-owl.net/uniprot.owl#"

def r = { String s ->
  if (s == "part-of") {
    factory.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/BFO_0000050"))
  } else if (s == "has-part") {
    factory.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/BFO_0000051"))
  } else {
    factory.getOWLObjectProperty(IRI.create("http://phenomebrowser.net/#"+s))
  }
}

def c = { String s ->
  factory.getOWLClass(IRI.create(s))
}

def addAnno = {resource, prop, cont ->
  OWLAnnotation anno = factory.getOWLAnnotation(
    factory.getOWLAnnotationProperty(prop.getIRI()),
    factory.getOWLLiteral(cont))
  def axiom = factory.getOWLAnnotationAssertionAxiom(resource.getIRI(),
                                                     anno)
  manager.addAxiom(ont,axiom)
}

def genericProteinNames = [:] // maps a name to an OWLClass
def proteinCounter = 0
new File("ttl").eachFile { ttlfile ->
  Model rdfModel = ModelFactory.createDefaultModel()
  rdfModel = rdfModel.read(new FileReader(ttlfile), "http://foobar#", "TURTLE")

  def p
  p = rdfModel.createProperty(up+"mnemonic")
  def label
  rdfModel.listStatements(null, p, null).each { stmt ->
    label = stmt.getObject()?.toString()
  }
  def iri
  rdfModel.listStatements(null, RDF.type, rdfModel.createResource(up+"Protein")).each { stmt ->
    iri = stmt.getSubject()?.toString()
  }

  def cl = c(iri) // create class
  addAnno(cl, OWLRDFVocabulary.RDFS_LABEL, label) // add label
  def cls = c(iri+"_all")
  addAnno(cls, OWLRDFVocabulary.RDFS_LABEL, "All $label in the universe") // add label
  manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(cl, fac.getOWLObjectSomeValuesFrom(r("member-of"), cls)))

  def genericNames = []
  def first = true
  def firstName
  rdfModel.listStatements(null, RDF.type, rdfModel.createResource(up+"Structured_Name")).each { stmt ->
    rdfModel.listStatements(stmt.getSubject(), rdfModel.createProperty(up+"fullName"), null).each { stmt2 ->
      if (first) {
	first = false
	firstName = stmt2.getObject()?.toString()
      }
      genericNames << stmt2.getObject()?.toString()
    }
  }
  if (! (firstName in genericProteinNames.keySet())) {
    def ncl = c(onturi+"GENERICPROTEIN_$proteinCounter")
    def ncla = c(onturi+"GENERICPROTEIN_$proteinCounter"+"_all")
    proteinCounter += 1
    genericNames.each {
      addAnno(ncl, OWLRDFVocabulary.RDFS_LABEL, it) // add label
      addAnno(ncla, OWLRDFVocabulary.RDFS_LABEL, "All $it in the universe") // add label
      genericProteinNames[it] = ncl
      genericProteinNames["All $it in the universe"] = ncla
    }
  }
  def ncl = genericProteinNames[firstName]
  def ncla = genericProteinNames["All $firstName in the universe"]
  manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(cl, ncl))
  manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(cls, ncla))

  def isoforms = new LinkedHashSet()
  rdfModel.listStatements(null, rdfModel.createProperty(up+"sequence"), null).each { stmt ->
    isoforms.add(stmt.getObject()?.toString())
  }
  isoforms.each { iso ->
    def icl = c(iso)
    addAnno(icl, OWLRDFVocabulary.RDFS_LABEL, "Isoform of $label") // add label
    manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(icl, fac.getOWLObjectSomeValuesFrom(r("isoform-of"), cl)))
    manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(icl, cl))
  }
  
  def gofunctions = new LinkedHashSet()
  rdfModel.listStatements(rdfModel.createResource(iri), rdfModel.createProperty(up+"classifiedWith"), null).each { stmt ->
    if(stmt.getObject()?.toString().startsWith("http://purl.obolibrary.org/obo/")) {
      gofunctions.add(stmt.getObject()?.toString())
    }
  }
  gofunctions.each { fun ->
    manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(cls, fac.getOWLObjectSomeValuesFrom(r("has-member"), fac.getOWLObjectSomeValuesFrom(r("has-function"), c(fun)))))
  }

  def organism
  rdfModel.listStatements(rdfModel.createResource(iri), rdfModel.createProperty(up+"organism"), null).each { stmt ->
    organism = stmt.getObject()?.toString()
  }
  manager.addAxiom(ont, factory.getOWLSubClassOfAxiom(cl, fac.getOWLObjectSomeValuesFrom(r("created-in-organism"), c(organism))))

  /*
  def tissues = new LinkedHashSet()
  rdfModel.listStatements(rdfModel.createResource(iri), rdfModel.createProperty(up+"isolatedFrom"), null).each { stmt ->
    tissues.add(stmt.getObject()?.toString())
  }
  println tissues
  */
}

manager.saveOntology(ont, IRI.create("file:"+args[0]))
