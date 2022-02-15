from paulssonlab.cloning.sequence import DsSeqRecord
from paulssonlab.cloning.workflow import normalize_seq_upper
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
import xmltodict
from lxml import etree
from lxml.builder import E
from cytoolz import get_in
from datetime import datetime
import time
import re

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
    ####
    # 'molecule_type'
    ####
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


def dumps_geneious(seq, sep="; "):
    root = dump_geneious(seq, sep=sep)
    return etree.tostring(seq, encoding="utf8")


def dump_geneious(seq, sep="; "):
    # SeqRecord/Seq
    if isinstance(seq, SeqRecord):
        # geneious interprets genbank dates (e.g., "29-OCT-2018") in local time zone
        # we instead correctly interpret them as UTC
        date_str = seq.annotations.get("date")
        if date_str:
            epoch = datetime(1970, 1, 1, 0, 0, 0)
            timestamp = str(
                int(
                    (datetime.strptime(date_str, "%d-%b-%Y") - epoch).total_seconds()
                    * 1000
                )
            )
        else:
            timestamp = str(int(time.time() * 1000))
        fields = []
        stored_fields = [
            E.standardField(code=name)
            for name in ("modified_date", "molType", "topology", "geneticCode")
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
            is_circular = "true"
        else:
            is_circular = "false"
        fields.append(E.isCircular(is_circular))
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
        # not sure a '.' is the right default value if missing genbank ID
        # that might just be bad handling by geneious or biopython (when parsing/outputting genbank)
        urn = E.urn(f"urn:sequence:genbank:{gid or '.'}", type="urn")
        created = E.created(timestamp, type="date")
        name = E.name(seq.name or "")
        description = E.description(seq.description or "")
        sequence_annotations = _dump_geneious_features(seq.features, sep=sep)
        char_sequence = E.charSequence(normalize_seq_upper(seq))
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
        xml = E.XMLSerialisableRootElement(
            E.fields(*fields),
            urn,
            created,
            E.storedFields(*stored_fields),
            name,
            description,
            sequence_annotations,
            char_sequence,
            E.INSD_originalElements(*insd_originalelements),
            # annotationsRevision="1",
            # charSequenceValidated="true",
        )
        return xml
    elif isinstance(seq, Seq):
        return dump_geneious(DsSeqRecord(seq))
    # primer
    elif isinstance(seq, Primer):
        pass
    else:
        raise ValueError(f"cannot dump type: {type(seq)}")


def _dump_geneious_features(features, sep=" / "):
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
                description = E.description(sep.join(value))
                continue  # geneious encodes the label qualifier in <description> instead
            qualifier = E.qualifier(E.name(name), E.value(sep.join(value)))
            qualifiers.append(qualifier)
        annotation_elements = [
            E.type(feature.type),
            E.intervals(*intervals),
            E.qualifiers(*qualifiers),
        ]
        if description is not None:
            # if label qualifier is missing, geneious seems to use f"{type} {qualifiers['source']}"
            # as the description. we don't do that
            annotation_elements = [description, *annotation_elements]
        element = E.annotation(*annotation_elements)
        elements.append(element)
    return E.sequenceAnnotations(*elements)


class GeneiousClient:
    # HOW DO WE PLAN TO INTEGRATE THIS WITH GDRIVE SYNCING?
    # what are the different document types (children of Folder)?
    # can geneious handle DsSeqRecord overhangs?
    # RENAME clear_cache -> rollback, save -> commit
    # refactor GDriveClient into a mixin

    # .raw for {"mimeType": ..., "content": hhh, "xml": "..."}
    # .xml for xml's
    # [] or .content for SeqRecords (how to handle duplicate names?? only error on getitem)

    def __init__(self, conn):
        self.conn = conn
        self.rollback()

    def rollback(self):
        pass

    @classmethod
    def connect(self, **kwargs):
        db_url = sqlalchemy.engine.URL.create(**kwargs)
        engine = sqlalchemy.create_engine(db_url)
        conn = engine.connect()
        self(conn)

    def __contains__(self, key):
        pass

    def __iter__(self):
        return iter(self.keys())

    def keys(self):
        pass

    def items(self):
        for key in self:
            yield key, self[key]

    def __getitem__(self, key):
        pass

    def __setitem__(self, key, value):
        pass

    def __delitem__(self, key):
        pass

    def _raw_contains(self, key):
        pass

    def _raw_delitem(self, key):
        pass

    def _raw_iter(self, key):
        pass

    def update(self, other):
        for key, value in other.items():
            self[key] = value

    def commit(self):
        pass
        # start transaction
        # read next_table_id
        # make new folders
        # delete empty folders
        # add/delete documents
        # update next_table_id
        # commit transaction, retry if failed?
        ################################################################################
        keys_to_trash = set()
        # REMOVE files/folders = None in local
        for key, value in self.local.items():
            if value is None:
                if key in self.remote or key in self.remote_folders:
                    keys_to_trash.add(key)
        for key in keys_to_trash:
            file = self.remote.get(key)
            is_folder = False
            if file is None:
                file = self.remote_folders.get(key)
                is_folder = True
            if file is None:
                continue
            self.client.files().update(
                fileId=file["id"],
                addParents=self._trash_folder_id(),
                removeParents=",".join(file["parents"]),
            ).execute()
            if is_folder:
                # remove all remote files/folders with a key prefixed by "key"
                for k in self.remote:
                    if key == k[: len(key)]:
                        del self.remote[k]
                for k in self.remote_folders:
                    if key == k[: len(key)]:
                        del self.remote_folders[k]
            else:
                del self.remote[key]
        # MAKE FOLDERS for all parent keys for all local keys
        # list() makes a copy so we can delete keys as we iterate
        for key, value in list(self.local.items()):
            if value is not None:
                for i in range(1, len(key)):
                    folder_key = key[:i]
                    if folder_key not in self.remote_folders:
                        folder = make_drive_folder(
                            self.client,
                            folder_key[-1],
                            self.remote_folders[folder_key[:-1]]["id"],
                        )
                        self.remote_folders[folder_key] = folder
                if key in self.remote:
                    if not overwrite:
                        raise ValueError(f"file already exists remotely: {key}")
                    file_id = self.remote[key]["id"]
                else:
                    file_id = None
                upload_drive(
                    self.client,
                    self.bytes[key],
                    key[-1],
                    mimetype=value["mimeType"],
                    file_id=file_id,
                    parent=self.remote_folders[key[:-1]]["id"],
                )
                self.remote[key] = self.local[key]
                del self.local[key]
        # TRASH ALL remaining folders
        if remove_empty_folders:
            nonempty_folders = set([()])  # never remove root
            for key in self.remote:
                for i in range(1, len(key)):
                    nonempty_folders.add(key[:i])
            empty_folders = set(self.remote_folders) - nonempty_folders
            for key in empty_folders:
                folder = self.remote_folders[key]
                self.client.files().update(
                    fileId=folder["id"],
                    addParents=self._trash_folder_id(),
                    removeParents=",".join(folder["parents"]),
                ).execute()
                del self.remote_folders[key]
