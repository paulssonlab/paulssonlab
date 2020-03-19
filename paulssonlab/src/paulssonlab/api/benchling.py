import benchlingapi


def get_folder(session, path=None, project=None, project_id=None):
    if project and project_id:
        raise ValueError("cannot specify both project and project_id")
    if project:
        if not isinstance(project, benchlingapi.models.Project):
            project = session.Project.find_by_name(project)
        project_id = project.id
    top_level_folders = session.Folder.list(
        parent_folder_id="NO_PARENT", project_id=project_id
    )
    if len(top_level_folders) != 1:
        raise Exception("expecting only one top level folder")
    folder = top_level_folders[0]
    if path:
        path_components = path.split("/")
        parent_folder_id = "NO_PARENT"
        for p in path_components:
            if not p:
                continue
            parent_folder_id = folder.id
            folder = session.Folder.find_by_name(
                p, parent_folder_id=parent_folder_id, project_id=project_id
            )
    return folder


def genbank_to_benchling(
    session,
    gb,
    name,
    folder_id,
    accession=None,
    long_annotations=True,
    custom_fields=None,
    sep=" / ",
):
    if gb.annotations["molecule_type"] != "ds-DNA":
        raise ValueError(
            f"unexpected value for molecule_type: {gb.annotations['molecule_type']}"
        )
    bases = gb.seq
    annotations = []
    translations = []
    for feature in gb.features:
        if feature.type == "source":
            feature_name = "source"
        else:
            if long_annotations:
                feature_name = sep.join(feature.qualifiers["label"])
                note = feature.qualifiers.get("note", None)
                if note is not None:
                    note = sep.join(note)
                if note:
                    feature_name += f" ({note})"
                gene = feature.qualifiers.get("gene", None)
                if gene is not None:
                    gene = sep.join(gene)
                if gene:
                    feature_name += f" (gene: {gene})"
                product = feature.qualifiers.get("product", None)
                if product is not None:
                    product = sep.join(product)
                if product:
                    feature_name += f" (product: {product})"

            else:
                feature_name = feature.qualifiers["label"] or feature.qualifiers["note"]
        if len(feature.location.parts) > 1:
            feature_name_all = feature_name + " [all]"
        else:
            feature_name_all = feature_name
        if len(feature.location.parts) > 1:
            for loc in feature.location.parts:
                start = int(loc.start)
                end = int(loc.end)
                annotation = {
                    "start": start,
                    "end": end,
                    "strand": loc.strand,
                    "name": feature_name + f" [{start}-{end}]",
                    "type": feature.type,
                }
                annotations.append(annotation)
        annotation = {
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "strand": feature.location.strand,
            "name": feature_name_all,
            "type": feature.type,
        }
        annotations.append(annotation)
        if feature.type == "CDS":  # TODO: are there any more CDS-like types?
            if int(feature.qualifiers["codon_start"][0]) != 1:
                raise ValueError("cannot handle codon_start != 1")
            translation = {
                "start": annotation["start"],
                "end": annotation["end"],
                "strand": annotation["strand"],
                "aminoAcids": feature.qualifiers["translation"],
                "regions": [
                    {"start": int(loc.start), "end": int(loc.end)}
                    for loc in feature.location.parts
                ],
            }
            translations.append(translation)
    translations = []
    if gb.annotations["topology"] == "circular":
        is_circular = True
    elif gb.annotations["topology"] == "linear":
        is_circular = False
    else:
        raise ValueError(f"unexpected value for topology: {gb.annotations['topology']}")
    _custom_fields = {}
    _custom_fields["organism"] = gb.annotations["organism"]
    _custom_fields["source"] = gb.annotations["source"]
    _custom_fields["division"] = gb.annotations["data_file_division"]
    _custom_fields["keywords"] = ",".join(gb.annotations["keywords"])
    _custom_fields["definition"] = gb.description
    if accession is None:
        accession = ",".join(gb.annotations["accessions"])
        if accession and accession != ".":
            _custom_fields["accession"] = accession
    else:
        _custom_fields["accession"] = accession
    if custom_fields:
        _custom_fields = {**_custom_fields, **custom_fields}
    _custom_fields = {k: {"value": v} for k, v in _custom_fields.items()}
    dna = session.DNASequence(
        name=name,
        folder_id=folder_id,
        bases=bases,
        annotations=annotations,
        translations=translations,
        is_circular=is_circular,
        custom_fields=_custom_fields,
    )
    return dna
