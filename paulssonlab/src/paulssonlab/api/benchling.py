import benchlingapi
from paulssonlab.cloning.sequence import get_seq
import json
import hashlib
import jsondiff
from copy import deepcopy

HASH_FIELD = "_hash"

# TODO: horribly hacky monkey-patch
benchlingapi.models.Folder.CREATE_SCHEMA = dict(
    only=("id", "name", "parent_folder_id", "project_id", "archive_record")
)
benchlingapi.models.Oligo.UPDATE_SCHEMA = benchlingapi.models.Oligo.CREATE_SCHEMA
# handle autoreloading
if benchlingapi.models.mixins.CreateMixin not in benchlingapi.models.Folder.__bases__:
    benchlingapi.models.Folder.__bases__ = (
        benchlingapi.models.mixins.CreateMixin,
    ) + benchlingapi.models.Folder.__bases__
if benchlingapi.models.mixins.ListMixin not in benchlingapi.models.Oligo.__bases__:
    benchlingapi.models.Oligo.__bases__ = (
        benchlingapi.models.mixins.ListMixin,
    ) + benchlingapi.models.Oligo.__bases__
if benchlingapi.models.mixins.UpdateMixin not in benchlingapi.models.Oligo.__bases__:
    benchlingapi.models.Oligo.__bases__ = (
        benchlingapi.models.mixins.UpdateMixin,
    ) + benchlingapi.models.Oligo.__bases__


def get_project_root(session, project_name):
    project = session.Project.find_by_name(project_name)
    folders = session.Folder.list(project_id=project.id, parent_folder_id="NO_PARENT")
    root_folder = folders[0]
    return root_folder


def _dictsorter(x):
    return tuple(x[k] for k in sorted(x.keys()))


def strip_dnasequence(seq):
    seq = deepcopy(seq)
    for key in ("id", "creator", "webUrl"):
        seq.pop(key, None)
    if "translations" in seq:
        for translation in seq["translations"]:
            for key in ("start", "end"):
                translation.pop(key, None)
    if "annotations" in seq:
        for annotation in seq["annotations"]:
            annotation.pop("color", None)
        seq["annotations"] = sorted(seq["annotations"], key=_dictsorter)
    if "customFields" in seq:
        seq["customFields"].pop(HASH_FIELD, None)
    return seq


def diff_dnasequence(a, b):
    for d in (a, b):
        strip_dnasequence(d)
    diff = jsondiff.diff(a, b)
    return diff


def hash_dnasequence(seq):
    seq = strip_dnasequence(seq)
    s = json.dumps(seq, sort_keys=True)
    # TODO: add a usedforsecurity=False flag to the constructor
    # in Python 3.9
    return hashlib.sha256(s.encode()).hexdigest()


def upload_sequence(
    root_folder,
    path,
    seq,
    oligo=False,
    overwrite=True,
    check_existing=True,
    update_in_place=True,
    use_hash=True,
    cache=None,
):
    session = root_folder.session
    if cache is None:
        cache = {}
    if isinstance(path, str):
        path = (path,)
    if len(path) >= 2:
        root_folder = ensure_folder(root_folder, path[:-1], cache=cache)
    will_update = False
    if check_existing:
        if oligo:
            entity_cls = session.Oligo
        else:
            entity_cls = session.DNASequence
        existing_dna = entity_cls.find_by_name(path[-1], folder_id=root_folder.id)
        if existing_dna is not None:
            if overwrite is True:
                if update_in_place:
                    will_update = True
                else:
                    existing_dna.archive()
            elif overwrite is False:
                raise ValueError(f"Benchling sequence already exists: {path}")
            elif overwrite is None:
                return existing_dna
    if oligo:
        bases = str(get_seq(seq))
        # we set a value for lengths and add empty aliases, fields, custom_fields
        # so that we can compute diffs against server-side models
        dna = session.Oligo(
            name=path[-1],
            bases=bases,
            length=len(bases),
            fields={},
            aliases=[],
            custom_fields={},
            folder_id=root_folder.id,
        )
    else:
        dna = seqrecord_to_benchling(session, root_folder.id, path[-1], seq)
    if will_update:
        dna_data = dna.dump()
        existing_dna_data = existing_dna.dump()
        if not oligo:
            new_hash = dna.custom_fields[HASH_FIELD]["value"]
            try:
                old_hash = existing_dna.custom_fields[HASH_FIELD]["value"]
            except:
                old_hash = None
            if use_hash and old_hash == new_hash:
                # don't update anything
                return existing_dna
        if len(diff_dnasequence(dna_data, existing_dna_data)) != 0:
            dna.id = existing_dna.id
            dna.update()
    else:
        dna.save()
    return dna


def ensure_folder(root_folder, path, cache=None):
    session = root_folder.session
    if cache is None:
        cache = {}
    if isinstance(path, str):
        path = (path,)
    for idx in range(1, len(path) + 1):
        current_path = path[:idx]
        if current_path in cache:
            root_folder = cache[current_path]
        else:
            new_root_folder = session.Folder.find_by_name(
                current_path[-1], parent_folder_id=root_folder.id
            )
            if new_root_folder is None:
                new_root_folder = session.Folder(
                    name=current_path[-1], parent_folder_id=root_folder.id
                )
                new_root_folder.save()
            root_folder = cache[current_path] = new_root_folder
    return root_folder


def seqrecord_to_benchling(
    session, folder_id, name, seq, long_annotations=True, custom_fields=None, sep=" / "
):
    molecule_type = seq.annotations.get("molecule_type")
    if molecule_type and molecule_type != "ds-DNA":
        raise ValueError(
            f"unexpected value for molecule_type: {seq.annotations['molecule_type']}"
        )
    bases = seq.seq
    length = len(bases)
    annotations = []
    translations = []
    for feature in seq.features:
        if feature.type == "source":
            feature_name = "source"
        else:
            if long_annotations:
                feature_names = []
                label = feature.qualifiers.get("label")
                if label is not None:
                    label = sep.join(label)
                    feature_names.append(label)
                note = feature.qualifiers.get("note")
                if note is not None:
                    note = sep.join(note)
                    feature_names.append(note)
                gene = feature.qualifiers.get("gene")
                if gene is not None:
                    gene = sep.join(gene)
                    feature_names.append(f"gene: {gene}")
                product = feature.qualifiers.get("product")
                if product is not None:
                    product = sep.join(product)
                    feature_names.append(f"product: {product}")
                if not len(feature_names):
                    feature_name = ""
                else:
                    feature_name = " ".join(
                        [feature_names[0], *[f"({n})" for n in feature_names[1:]]]
                    )
            else:
                feature_name = (
                    feature.qualifiers["label"] or feature.qualifiers["note"] or ""
                )
        if len(feature.location.parts) > 1:
            feature_name_all = feature_name + " [all]"
        else:
            feature_name_all = feature_name
        if len(feature.location.parts) > 1:
            for loc in feature.location.parts:
                start = int(loc.start)
                end = int(loc.end)
                # modulo length here and below because benchling seems to encode
                # a whole-plasmid annotation as 0-0 instead of 0-end
                annotation = {
                    "start": start % length,
                    "end": end % length,
                    "strand": loc.strand,
                    "name": feature_name + f" [{start}-{end}]",
                    "type": feature.type,
                }
                annotations.append(annotation)
        annotation = {
            "start": int(feature.location.start) % length,
            "end": int(feature.location.end) % length,
            "strand": feature.location.strand,
            "name": feature_name_all,
            "type": feature.type,
        }
        annotations.append(annotation)
        if feature.type == "CDS":  # TODO: are there any more CDS-like types?
            if int(feature.qualifiers["codon_start"][0]) != 1:
                raise ValueError("cannot handle codon_start != 1")
            translation = {
                "strand": annotation["strand"],
                "aminoAcids": feature.qualifiers["translation"],
                "regions": [
                    {"start": int(loc.start) % length, "end": int(loc.end) % length}
                    for loc in feature.location.parts
                ],
            }
            translations.append(translation)
    if seq.annotations["topology"] == "circular":
        is_circular = True
    elif seq.annotations["topology"] == "linear":
        is_circular = False
    else:
        raise ValueError(
            f"unexpected value for topology: {seq.annotations['topology']}"
        )
    _custom_fields = {}
    if "organism" in seq.annotations:
        _custom_fields["organism"] = seq.annotations["organism"]
    if "source" in seq.annotations:
        _custom_fields["source"] = seq.annotations["source"]
    if "data_file_division" in seq.annotations:
        _custom_fields["division"] = seq.annotations["data_file_division"]
    if "keywords" in seq.annotations:
        _custom_fields["keywords"] = ",".join(seq.annotations["keywords"])
    _custom_fields["definition"] = seq.description
    if "accessions" in seq.annotations:
        accession = ",".join(seq.annotations["accessions"])
        if accession and accession != ".":
            _custom_fields["accession"] = accession
    if custom_fields:
        _custom_fields = {**_custom_fields, **custom_fields}
    _custom_fields = {k: {"value": v} for k, v in _custom_fields.items()}
    # we set a value for lengths and add empty aliases, fields, primers
    # so that we can compute diffs against server-side models
    dna = session.DNASequence(
        name=name,
        folder_id=folder_id,
        bases=bases,
        length=length,
        aliases=[],
        fields={},
        primers=[],
        annotations=annotations,
        translations=translations,
        is_circular=is_circular,
        custom_fields=_custom_fields,
    )
    dna.custom_fields["_hash"] = {"value": hash_dnasequence(dna.dump())}
    return dna
