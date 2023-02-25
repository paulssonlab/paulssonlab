from paulssonlab.cloning.sequence import DsSeqRecord, normalize_seq_upper
from paulssonlab.cloning.primers import Primer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    CompoundLocation,
    ExactPosition,
    BeforePosition,
    AfterPosition,
    Reference,
)
from Bio import GenBank
import sqlalchemy
from lxml import etree
from lxml.builder import E
from cytoolz import get_in
from datetime import datetime
import time
import re
import uuid
from copy import deepcopy

GENEIOUS_TO_GENBANK = {
    "fields/topology": "topology",
    "INSD_originalElements/INSDSeq_division": "data_file_division",
    "fields/gid": "gi",
    "INSD_originalElements/INSDSeq_source": "source",
    "fields/organism": "organism",
}

INSD_TO_REFERENCE = {
    "INSDReference_consortium": "consrtm",
    "INSDReference_title": "title",
    "INSDReference_journal": "journal",
    "INSDReference_pubmed": "pubmed_id",
}


def make_urn(location="local", user="paulssonlab"):
    s = uuid.uuid4().hex
    return f"urn:{location}:{user}:{s[:2]}-{s[2:9]}"


def geneious_timestamp():
    return int(time.time() * 1000)


def connect(**kwargs):
    db_url = sqlalchemy.engine.URL.create(**kwargs)
    engine = sqlalchemy.create_engine(db_url, future=True)
    sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine, future=True)
    return sessionmaker


def loads_geneious(s):
    return load_geneious(etree.fromstring(s))


def _load_geneious_reference(xml, seq_length):
    reference = Reference()
    for src, dest in INSD_TO_REFERENCE.items():
        element = xml.find(src)
        if element is not None:
            setattr(reference, dest, element.text)
    authors = [a.text for a in xml.findall("INSDReference_authors/INSDAuthor")]
    if authors:
        reference.authors = authors
    position_element = xml.find("INSDReference_position")
    if position_element is not None:
        # we do not handle compound locations, etc.
        # (BioPython's hacky GenBank parser makes it hard to do this)
        reference.location = GenBank._loc(position_element.text, seq_length, 0)
    return reference


def load_geneious(xml, sep="; "):
    fields = xml.find("fields")
    override_document_type = fields.find("overrideDocumentType")
    if (
        override_document_type
        and override_document_type.text
        == "com.biomatters.geneious.publicapi.documents.sequence.DefaultSequenceListDocument$NucleotideSequenceList"
    ):
        # sequence list
        seqs = [_load_geneious(d) for d in xml.findall("nucleotideSequence")]
        return seqs
    # single seq
    seq = Seq(xml.find("charSequence").text)
    annotations = {}
    # fields
    name = xml.find("name").text
    description = xml.find("description").text
    fields = xml.find("fields")
    accession_element = fields.find("accession")
    if accession_element is not None:
        accession = accession_element.text
        if accession.count(".") == 1:
            accession, version = accession.split(".")
            version = int(version)
            annotations["sequence_version"] = version
            id_ = f"{accession}.{version}"
        else:
            id_ = accession
        annotations["accessions"] = [accession]
    else:
        id_ = None
    # features
    features = [
        _load_geneious_feature(a) for a in xml.findall("sequenceAnnotations/annotation")
    ]
    source_feature = None
    for feature in features:
        if feature.type == "source":
            source_feature = feature
            break
    # annotations
    molecule_type = "ds-DNA"
    moltype_element = xml.find("fields/molType")
    if moltype_element is not None:
        molecule_type = moltype_element.text
        strandedness_element = xml.find("INSD_originalElements/INSDSeq_strandedness")
        if "-" not in molecule_type and strandedness_element is not None:
            strandedness = strandedness_element.text
            if strandedness == "double":
                molecule_type = "ds-" + molecule_type
            elif strandedness == "single":
                molecule_type = "ss-" + molecule_type
            if strandedness == "mixed":
                molecule_type = "ms-" + molecule_type
    annotations["molecule_type"] = molecule_type
    circular = None
    is_circular_element = xml.find("fields/isCircular")
    if is_circular_element is not None:
        if is_circular_element.text == "true":
            circular = True
    for src, dest in GENEIOUS_TO_GENBANK.items():
        element = xml.find(src)
        if element is not None:
            annotations[dest] = element.text
    date_element = xml.find("fields/modified_date")
    if date_element is not None:
        timestamp = int(date_element.text) / 1000
        modified_date = datetime.fromtimestamp(timestamp)
        annotations["date"] = modified_date.strftime("%d-%b-%Y")
    taxonomy_element = xml.find("fields/taxonomy")
    if taxonomy_element is not None:
        taxonomy = re.split(r"\s*;?\s+", taxonomy_element.text)
        if taxonomy:
            annotations["taxonomy"] = taxonomy
    references = [
        _load_geneious_reference(ref, len(seq))
        for ref in xml.findall("INSD_originalElements/INSDSeq_references/INSDReference")
    ]
    if references:
        annotations["references"] = references
    keywords = [
        k.text
        for k in xml.findall("INSD_originalElements/INSDSeq_keywords/INSDSeq_keyword")
    ]
    if keywords:
        annotations["keywords"] = keywords
    # TODO: ds-DNA
    res = DsSeqRecord(
        seq,
        id=id_,
        name=name,
        description=description,
        features=features,
        annotations=annotations,
        circular=circular,
    )
    return res


def _load_geneious_feature(xml):
    parts = []
    for interval in xml.findall("intervals/interval"):
        start = int(interval.find("minimumIndex").text)
        stop = int(interval.find("maximumIndex").text)
        direction_element = interval.find("direction")
        if direction_element is not None:
            direction = direction_element.text
        else:
            direction = None
        if direction == "leftToRight":
            strand = 1
        elif direction == "rightToLeft":
            strand = -1
        else:
            strand = 0
        if interval.attrib.get("beginsBeforeMinimumIndex") == "true":
            start = BeforePosition(start)
        else:
            start = ExactPosition(start)
        if interval.attrib.get("endsAfterMaximumIndex") == "true":
            stop = AfterPosition(stop)
        else:
            stop = ExactPosition(stop)
        parts.append(FeatureLocation(start, stop, strand=strand))
    if len(parts) == 0:
        raise ValueError("feature missing location")
    elif len(parts) == 1:
        location = parts[0]
    else:
        location = CompoundLocation(*parts)
    type_ = xml.find("type").text
    # skip description, geneious sets description using the contents
    # of other qualifiers
    feature = SeqFeature(location, type=type_)
    for qualifier in xml.findall("qualifiers/qualifier"):
        name = qualifier.find("name").text
        value = qualifier.find("value").text
        if name == "NCBI Feature Key":
            # skip, this is already included in type_
            continue
        feature.qualifiers[name] = [value]
    return feature


def dumps_geneious(*args, pretty_print=False, **kwargs):
    xmls = dump_geneious(*args, **kwargs)
    return [
        etree.tostring(xml, encoding=str, pretty_print=pretty_print) for xml in xmls
    ]


def dump_geneious(
    seq,
    urn,
    document_type="genbank",
    name=None,
    description=None,
    timestamp=None,
    sep="; ",
):
    if document_type == "genbank":
        # SeqRecord/Seq
        seq = DsSeqRecord(seq)
        if name is None:
            name = seq.name or ""
        if description is None:
            description = seq.description or ""
        # geneious interprets genbank dates (e.g., "29-OCT-2018") in local time zone
        # we instead correctly interpret them as UTC
        if timestamp is None:
            date_str = seq.annotations.get("date")
            if date_str:
                epoch = datetime(1970, 1, 1, 0, 0, 0)
                timestamp = int(
                    (datetime.strptime(date_str, "%d-%b-%Y") - epoch).total_seconds()
                    * 1000
                )
            else:
                timestamp = geneious_timestamp()
        timestamp = str(timestamp)
        fields = []
        stored_fields = [
            E.standardField(code=code)
            for code in ("modified_date", "molType", "topology", "geneticCode")
        ]
        insd_originalelements = []
        source_feature = {}
        for feature in seq.features:
            if feature.type == "source":
                source_feature = feature
                break
        # TODO: geneious imports commonName, db_xref, and chromosome from source_feature
        # we decide to skip this
        gid = seq.annotations.get("gi")
        if gid:
            fields.append(E.gid(gid))
        organism = seq.annotations.get("organism")
        if organism is not None:
            fields.append(E.organism(organism))
            stored_fields.append(E.standardField(code="organism"))
        fields.append(E.hasBeenModified("false", type="boolean"))
        if seq.annotations.get("topology") == "circular":
            fields.append(E.isCircular("true"))
            # isCircular element should not exist for linear sequences!
            # if topology == "linear" but we add <isCircular>false</isCircular>
            # geneious will still show circular sequence viewer
        topology = seq.annotations.get("topology", "linear")
        if topology not in ("circular", "linear"):
            raise ValueError(f"do not understand topology '{topology}'")
        fields.append(E.topology(topology))
        fields.append(E.geneticCode("Standard"))
        # TODO: not sure what separator should be used for multiple accessions
        accessions = seq.annotations.get("accessions", [])
        if accessions:
            sequence_version = seq.annotations.get("sequence_version")
            if sequence_version:
                accession_format = "{}." + str(sequence_version)
            else:
                accession_format = "{}"
            # SEE: Bio.GenBank._split_accessions
            fields.append(
                E.accession(" ".join(accession_format.format(a) for a in accessions))
            )
            stored_fields.append(E.standardField(code="accession"))
        taxonomy = seq.annotations.get("taxonomy", [])
        if taxonomy:
            fields.append(E.taxonomy(sep.join(taxonomy)))
            stored_fields.append(E.standardField(code="taxonomy"))
        fields.append(E.modified_date(timestamp, type="date"))
        molecule_type = seq.annotations.get("molecule_type", "ds-DNA")
        if molecule_type.startswith("ss-"):
            strandedness = "single"
            molecule_type = molecule_type[3:]
        elif molecule_type.startswith("ms-"):
            strandedness = "mixed"
            molecule_type = molecule_type[3:]
        elif molecule_type.startswith("ds-"):
            strandedness = "double"
            molecule_type = molecule_type[3:]
        else:
            strandedness = None
        fields.append(E.molType(molecule_type))
        plugin_document_urn = f"urn:sequence:genbank:{gid or 'null'}"
        urn_element = E.urn(plugin_document_urn, type="urn")
        created_element = E.created(timestamp, type="date")
        name_element = E.name(name)
        description_element = E.description(description)
        sequence_annotations = _dump_geneious_features(seq.features, sep=sep)
        normalized_seq = normalize_seq_upper(seq)
        char_sequence = E.charSequence(normalized_seq)
        insd_originalelements.append(E.INSDSeq_length(str(len(seq))))
        if seq.name:
            insd_originalelements.append(E.INSDSeq_locus(seq.name))
        if strandedness:
            insd_originalelements.append(E.INSDSeq_strandedness(strandedness))
        division = seq.annotations.get("data_file_division")
        if division:
            insd_originalelements.append(E.INSDSeq_division(division))
        if gid:
            other_seqids = [E.INSDSeqid(f"gi|{gid}")]
        else:
            other_seqids = []
        insd_originalelements.append(getattr(E, "INSDSeq_other-seqids")(*other_seqids))
        keywords = [k for k in seq.annotations.get("keywords", []) if k]
        if keywords:
            insd_originalelements.append(
                E.INSDSeq_keywords(*[E.INSDSeq_keyword(k) for k in keywords])
            )
        source = seq.annotations.get("source")
        if not source and source_feature:
            # if organism is missing/empty, is there any other qualifier geneious pulls from?
            source = sep.join(feature.qualifiers.get("organism"))
        if source:
            insd_originalelements.append(E.INSDSeq_source(source))
        insd_references = []
        for idx, reference in enumerate(seq.annotations.get("references", [])):
            insd_reference = [E.INSDReference_reference(str(idx + 1))]
            if reference.location:
                loc = reference.location[0]
                # for more on indexing, see comment in _dump_geneious_features
                insd_reference.append(
                    E.INSDReference_position(f"{int(loc.start)+1}..{int(loc.end)}")
                )
            if reference.authors:
                authors = re.split(r",\s+|,?\s+and\s+", reference.authors)
                insd_reference.append(
                    E.INSDReference_authors(*[E.INSDAuthor(a) for a in authors])
                )
            if reference.consrtm:
                insd_reference.append(E.INSDReference_consortium(reference.consrtm))
            if reference.title:
                insd_reference.append(E.INSDReference_title(reference.title))
            if reference.journal:
                insd_reference.append(E.INSDReference_journal(reference.journal))
            if reference.pubmed_id:
                insd_reference.append(E.INSDReference_pubmed(reference.pubmed_id))
            insd_references.append(E.INSDReference(*insd_reference))
        insd_originalelements.append(E.INSDSeq_references(*insd_references))
        plugin_document_xml = E.XMLSerialisableRootElement(
            E.fields(*fields),
            urn_element,
            created_element,
            E.storedFields(*stored_fields),
            name_element,
            description_element,
            sequence_annotations,
            char_sequence,
            E.INSD_originalElements(*insd_originalelements),
            # annotationsRevision="1",
            # charSequenceValidated="true",
        )
        document_xml = dump_geneious_document_xml(
            plugin_document_xml=plugin_document_xml,
            plugin_document_urn=plugin_document_urn,
            urn=urn,
            name=name,
            description=description,
            genetic_code="Standard",
            molecule_type=molecule_type,
            topology=topology,
            organism=organism,
            seq=normalized_seq,
            timestamp=timestamp,
            document_class="com.biomatters.plugins.ncbi.documents.GenBankNucleotideSequence",
        )
        return document_xml, plugin_document_xml
    # primer
    elif document_type == "primer":
        if not isinstance(seq, Primer):
            seq = Primer(seq)
        if name is None:
            raise ValueError("name must be provided")
        if description is None:
            raise ValueError("description must be provided")
        if timestamp is None:
            raise ValueError("timestamp must be provided")
        timestamp = str(timestamp)
        normalized_seq = normalize_seq_upper(seq.seq)
        binding = normalize_seq_upper(seq.binding)
        overhang = normalize_seq_upper(seq.overhang)
        minimum_index = len(overhang) + 1
        maximum_index = len(seq)
        fields_element = E.fields(E.oligoType("Primer"), E.molType("DNA"))
        created_element = E.created(timestamp, type="date")
        stored_fields = [
            E.standardField(code="organism"),
            E.field(
                E.DocumentField(
                    E.name("Oligo Type"),
                    E.description(
                        "The type of oligonulceotide, one of: Primer and Probe"
                    ),
                    E.code("oligoType"),
                    E.type("java.lang.String"),
                    E.enumerationOptions(E.value("Primer"), E.value("Probe")),
                ),
                **{
                    "class": "com.biomatters.geneious.publicapi.documents.DocumentField"
                },
            ),
        ]
        annotation_uuid = str(uuid.uuid4())
        annotations = [
            E.annotation(
                E.description("Binding Region"),
                E.type("primer_bind"),
                E.intervals(
                    E.interval(
                        E.minimumIndex(str(minimum_index)),
                        E.maximumIndex(str(maximum_index)),
                        E.direction("leftToRight"),
                    )
                ),
                E.qualifiers(
                    E.qualifier(E.name("annotation group"), E.value(annotation_uuid)),
                    E.qualifier(E.name("Sequence"), E.value(binding)),
                    E.qualifier(E.name("Extension"), E.value(overhang)),
                ),
            )
        ]
        # not sure why sequenceAnnotations and primerAnnotations contain duplicated information
        # need to deepcopy annotations, otherwise it only appears in primerAnnotation
        plugin_document_xml = E.XMLSerialisableRootElement(
            fields_element,
            created_element,
            E.storedFields(*stored_fields),
            E.name(name),
            E.description(description),
            E.sequenceAnnotations(*deepcopy(annotations)),
            E.charSequence(normalized_seq),
            E.primerAnnotation(
                *annotations,
                **{
                    "class": "com.biomatters.geneious.publicapi.documents.sequence.SequenceAnnotation"
                },
            ),
        )
        document_xml = dump_geneious_document_xml(
            plugin_document_xml=plugin_document_xml,
            plugin_document_urn="",
            urn=urn,
            name=name,
            description=description,
            oligo_type="Primer",
            topology="linear",
            seq=normalized_seq,
            timestamp=timestamp,
            document_class="com.biomatters.geneious.publicapi.implementations.sequence.OligoSequenceDocument",
        )
        return document_xml, plugin_document_xml
    else:
        raise ValueError(f"cannot dump type: {document_type}")


def dump_geneious_document_xml(
    num_bytes=None,
    plugin_document_xml=None,
    plugin_document_urn="",
    urn="",
    name=None,
    description=None,
    oligo_type=None,
    genetic_code=None,
    molecule_type=None,
    topology=None,
    organism=None,
    seq=None,
    timestamp=None,
    document_class=None,
):
    if num_bytes is None:
        num_bytes = len(etree.tostring(plugin_document_xml, encoding="utf8").decode())
    hidden_fields = [
        E.cache_plugin_document_urn(plugin_document_urn),
        E.cache_urn(urn, type="urn"),
        E.nucleotideSequenceWithQualityCount("0", type="int"),
        E.cache_name(name),
        E.cache_created(timestamp, type="date"),
        E.trimmedSequencesCount("0", type="int"),
    ]
    if description is not None:
        hidden_fields.append(E.description(description))
    hidden_fields_element = E.hiddenFields(*hidden_fields)
    fields = [
        E.document_size(str(num_bytes), type="bytes"),
        E.sequence_length(str(len(seq)), type="int"),
        E.sequence_residues(seq[:500]),
        E.modified_date(timestamp, type="date"),
    ]
    if genetic_code is not None:
        fields.append(E.geneticCode(genetic_code))
    if molecule_type is not None:
        fields.append(E.molType(molecule_type))
    if topology is not None:
        fields.append(E.topology(topology))
    if organism is not None:
        fields.append(E.organism(organism))
    if oligo_type is not None:
        fields.append(E.oligoType(oligo_type))
    fields_element = E.fields(*fields)
    xml = E.document(
        hidden_fields_element,
        fields_element,
        version="1.3-11",
        revisionNumber="1",
        geneiousVersion="2022.0.2",
        geneiousVersionMinimum="7.1",
        PluginDocument_FormatLastChanged="7.1",
        PluginDocument_FormatLastExtended="8.1",
        PluginDocument_OldestVersionSerializableTo="6.0",
        **{"class": document_class},
    )
    return xml


def _dump_geneious_features(features, sep="; "):
    elements = []
    for feature in features:
        intervals = []
        for part in feature.location.parts:
            # we use start+1, end because geneious seems to use 1-indexing
            # and maximumIndex is inclusive
            # (so unlike Python indexing in both respects)
            interval_elements = [
                E.minimumIndex(str(int(part.start) + 1)),
                E.maximumIndex(str(int(part.end))),
            ]
            if part.strand == 1:
                direction = "leftToRight"
            elif part.strand == -1:
                direction = "rightToLeft"
            else:
                direction = None
            if direction:
                interval_elements.append(E.direction(direction))
            attributes = {}
            if isinstance(part.start, BeforePosition):
                attributes["beginsBeforeMinimumIndex"] = "true"
            if isinstance(part.end, AfterPosition):
                attributes["endsAfterMaximumIndex"] = "true"
            interval = E.interval(*interval_elements, **attributes)
            intervals.append(interval)
        qualifiers = []
        description = None
        for name, value in feature.qualifiers.items():
            if name == "label":
                description = sep.join(value)
                continue  # geneious encodes the label qualifier in <description> instead
            qualifier = E.qualifier(E.name(name), E.value(sep.join(value)))
            qualifiers.append(qualifier)
        if description is None:
            # mimic what Geneious does
            if feature.type == "source":
                description = feature.type
                label = sep.join(feature.qualifiers.get("organism", []))
                if label:
                    description += f" {label}"
            else:
                description = ""
                for qualifier in ("gene", "product", "note"):
                    label = sep.join(feature.qualifiers.get(qualifier, []))
                    if label:
                        description = f"{label} "
                        break
                description += feature.type
        element = E.annotation(
            E.description(description),
            E.type(feature.type),
            E.intervals(*intervals),
            E.qualifiers(*qualifiers),
        )
        elements.append(element)
    return E.sequenceAnnotations(*elements)
