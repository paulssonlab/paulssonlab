def ligate(seqs, method, linear=True):
    if method == "goldengate":
        # _check_seq_compatibility
        # _reverse_complement_overhangs
        # _5prime_overhang
        pass
    elif method == "gibson":
        pass
    else:
        raise ValueError("method should be one of: goldengate, gibson")
    if len(seqs) < 2:
        raise ValueError("need at least two sequences to assemble")
    seq1 = seqs[0]
    seq2 = seqs[1]
    seq1_rc = _reverse_complement_overhangs(seq1)
    seq2_rc = _reverse_complement_overhangs(seq2)
    if _check_seq_compatibility(seq1, seq2):
        pass
    elif _check_seq_compatibility(seq1, seq2_rc):
        seqs[1] = seq2_rc
    elif _check_seq_compatibility(seq1_rc, seq2):
        seqs[0] = seq1_rc
    elif _check_seq_compatibility(seq1_rc, seq2_rc):
        seqs[0] = seq1_rc
        seqs[1] = seq2_rc
    else:
        raise ValueError(f"overhang mismatch when assembling sequences 0 and 1")
    seqs = [*seqs, None]
    seqs_to_join = []
    for idx in range(len(seqs) - 1):
        seq1 = seqs[idx]
        seq2 = seqs[idx + 1]
        if seq2 is not None:
            seq2_rc = _reverse_complement_overhangs(seq2)
            if _check_seq_compatibility(seq1, seq2):
                pass
            elif _check_seq_compatibility(seq1, seq2_rc):
                seq2 = seqs[idx + 1] = seq2_rc
            else:
                raise ValueError(
                    f"overhang mismatch when assembling sequences {idx} and {idx + 1}: {seq1[2]} does not match {seq2[1]} or {seq2_rc[1]}"
                )
        seqs_to_join.append(_5prime_overhang(seq1[1]))
        seqs_to_join.append(seq1[0])
        if seq2 is None:
            if linear:
                seqs_to_join.append(_5prime_overhang(seq1[2]))
            else:
                if not _check_seq_compatibility(seq1, seqs[0]):
                    raise ValueError(
                        f"overhang mismatch when assembling sequences {idx} and 0 (circularizing): {seq1[2]} does not match {seqs[0][1]}"
                    )
    # copy SeqRecords, add annotations for each part?? (including overhangs)
    joined_seq = join_seqs(seqs_to_join)
    joined_seq.annotations = {
        "molecule_type": "ds-DNA",
        "data_file_division": "SYN",
        "date": datetime.now(timezone.utc).isoformat(),
        "source": "synthetic DNA construct",
        "organism": "synthetic DNA construct",
    }
    if linear is False:
        joined_seq.annotations["topology"] = "circular"
    elif linear:
        joined_seq.annotations["topology"] = "linear"
    # TODO: add circular annotation?
    return joined_seq


def assemble(parts, linear=True):
    seqs_to_assemble = []
    for part in parts:
        if len(part) < 2 or len(part) >= 5:
            raise ValueError(
                "expected between two and four arguments: sequence, enzyme, part_name, part_type"
            )
        part_seq, enzyme, part_name, part_type = [*part, *[None] * (4 - len(part))]
        subseqs = re_digest(part_seq, enzyme, linear=linear)
        part_seq, overhang1, overhang2 = subseqs[0]
        if hasattr(part_seq, "features"):
            part_seq = deepcopy(part_seq)
            if part_type is None:
                part_type = "misc_feature"
            if part_name is not None:
                label = SeqFeature(
                    FeatureLocation(
                        -len(overhang1[0]), len(part_seq) + len(overhang2[0])
                    ),
                    type=part_type,
                )
                label.qualifiers["label"] = [part_name]
                features = [
                    feature for feature in part.features if feature.type != "source"
                ]
                part_seq.features = [label, *features]
        seqs_to_assemble.append((part_seq, overhang1, overhang2))
    assembly = ligate(seqs_to_assemble, linear=linear)
    return assembly
