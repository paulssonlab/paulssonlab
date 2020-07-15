import Bio


def format_well_name(plate, row, column):
    return f"{plate if plate and plate != 1 else ''}{row}{column}"


def _move_feature(feature, start, end):
    return SeqFeature(
        FeatureLocation(
            ExactPosition(start), ExactPosition(end), strand=feature.location.strand
        ),
        type=feature.type,
        id=feature.id,
        qualifiers=feature.qualifiers,
    )


def slice_feature(feature, start, end):
    f_start = feature.location.nofuzzy_start
    f_end = feature.location.nofuzzy_end
    endpoints = sorted([f_start, f_end, start, end])
    # TODO: check equality
    if ((f_end - f_start) >= 0) == ((end - start) >= 0):
        es = (endpoints[1:3],)
    else:
        es = (endpoints[:2], endpoints[2:])
    return [_move_feature(feature, e1, e2) for e1, e2 in es if e1 < e2]


def slice_seq(seq, start, end):
    if start is None:
        start = 0
    if end is None:
        end = len(seq)
    # TODO: handle circular slices
    if end < start:
        slice1 = slice_seq(seq, start, None)
        slice2 = slice_seq(seq, None, end)
        new_seq = slice1 + slice2
        # copy letter annotations
    else:
        new_seq = seq[start:end]
        if hasattr(seq, "features"):
            features = []
            for feature in seq.features:
                if (
                    end <= feature.location.nofuzzy_start
                    or start >= feature.location.nofuzzy_end
                ):
                    continue
                new_feature = feature._shift(-start)
                start_loc = new_feature.location.start
                end_loc = new_feature.location.end
                if start > feature.location.nofuzzy_start:
                    start_loc = ExactPosition(0)
                if end < feature.location.nofuzzy_end:
                    end_loc = ExactPosition(end - start)
                new_feature.location = FeatureLocation(
                    start_loc, end_loc, strand=new_feature.location.strand
                )
                features.append(new_feature)
            new_seq.features = features
    return new_seq
