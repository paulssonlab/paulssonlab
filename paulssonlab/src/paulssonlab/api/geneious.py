from paulssonlab.cloning.sequence import DsSeqRecord
from paulssonlab.cloning.workflow import normalize_seq_upper
from paulssonlab.cloning.primers import Primer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import sqlalchemy
import xmltodict
from lxml import etree
from lxml.builder import E
from cytoolz import get_in
import time


# def _parse_xml(s):
#     # use native dicts, since they are ordered starting in Python 3.6
#     return xmltodict.parse(s, dict_constructor=dict)


def load_geneious(s):
    return _load_geneious(_parse_xml(s))


def _load_geneious(xml):
    if "XMLSerialisableRootElement" in xml:
        xml = xml["XMLSerialisableRootElement"]
    fields = xml["fields"]
    if (
        fields.get("overrideDocumentType")
        == "com.biomatters.geneious.publicapi.documents.sequence.DefaultSequenceListDocument$NucleotideSequenceList"
    ):
        # sequence list
        seqs = [_load_geneious(d) for d in xml["nucleotideSequence"]]
        return seqs
    else:
        # single seq
        seq = Seq(xml["charSequence"])
        # fields
        name = xml.get("name")
        description = xml.get("description")
        # features
        features = [
            _load_geneious_feature(a)
            for a in get_in(("sequenceAnnotations", "annotation"), d)
        ]
        # annotations
        insd_annotations = xml.get("INSD_originalElements")
        if insd_annotations is not None:
            # TODO: this should be name
            id_ = insd_annotations["INSDSeq_locus"]
            annotations = {}
            for key in ("division", "source"):
                insd_key = f"INSDSeq_{key}"
                if insd_key in insd_annotations:
                    annotations[key] = insd_annotations[insd_key]
            if "INSDSeq_keywords" in insd_annotations:
                # TODO: only handles a single keyword here... is there ever more than one?
                annotations["keywords"] = [
                    get_in(("INSDSeq_keywords", "INSDSeq_keyword"), insd_annotations)
                ]
            if "INSDSeq_division" in annotations:
                pass
            # INSDSeq_references
        else:
            id_ = None
            annotations = None
        # TODO: circularity
        # TODO: ds-DNA
        res = DsSeqRecord(
            seq,
            id=id_,
            name=name,
            description=description,
            features=features,
            annotations=annotations,
        )
        return res


def _load_geneious_feature(d):
    intervals = d["intervals"]
    if "interval" in intervals:
        interval = intervals["interval"]
        start = interval["minimumIndex"]
        stop = interval["maximumIndex"]
        direction = interval["direction"]
        # TODO: strand
    else:
        raise ValueError(f"expecting 'interval' in annotation location {intervals}")
    type_ = d.get("type")
    # description = d.get("description")
    feature = SeqFeature(FeatureLocation(start, stop), type=type_)
    for qualifier in d["qualifiers"]:
        if qualifier["name"] == "NCBI Feature Key":
            # skip, this is already included in type_
            continue
        feature.qualifiers[qualifier["name"]] = [qualifier["value"]]
    return feature


def dump_geneious(seq, sep=" / "):
    # TODO: seq.id, oligos
    # SeqRecord/Seq
    if isinstance(seq, SeqRecord):
        now = str(int(time.time()))
        fields = []
        stored_fields = [
            E.standardField(code=name)
            for name in ("modified_date", "molType", "topology", "geneticCode")
        ]
        insd_originalelements = []
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
        accessions = seq.annotations.get("accessions", [])
        if accessions:
            fields.append(E.accession(sep.join(accessions)))
            stored_fields.append(E.standardField(code="accession"))
        taxonomy = seq.annotations.get("taxonomy", [])
        if taxonomy:
            fields.append(E.taxonomy(sep.join(taxonomy)))
            stored_fields.append(E.standardField(code="taxonomy"))
        fields.append(E.modified_date(now, type="date"))
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
        urn = E.urn("urn:sequence:genbank:.", type="urn")
        created = E.created(now, type="date")
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
        insd_originalelements.append(getattr(E, "INSDSeq_other-seqids"))
        insd_originalelements.append(
            E.INSDSeq_keywords(
                *[E.INSDSeq_keyword(k) for k in seq.annotations.get("keywords", [])]
            )
        )
        source = None
        for feature in seq.features:
            if feature.type == "source":
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
                insd_reference.append(
                    E.INSDReference_authors(E.INSDAuthor(reference.authors))
                )
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
            if part.strand == 1:
                direction = "leftToRight"
            elif part.strand == -1:
                direction = "rightToLeft"
            else:
                direction = "none"
            # we use start+1, end because geneious seems to use 1-indexing
            # and maximumIndex is inclusive
            # (so unlike Python indexing in both respects)
            interval = E.interval(
                E.minimumIndex(str(int(part.start) + 1)),
                E.maximumIndex(str(int(part.end))),
                E.direction(direction),
            )
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
