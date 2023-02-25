import dnachisel as dc
from Bio.Seq import Seq


def dnachisel_constraints_for_twist(
    seq,
    cds_location=None,
    avoid_enzymes=None,
    aa_to_codons=None,
    genetic_table="Bacterial",
    codon_optimization="use_best_codon",
):
    if aa_to_codons is None:
        aa_to_codons = codon.codons_by_relative_frequency()
    # TODO: could generate a genetic_table from codon_usage_table, if supplied
    constraints = []
    objectives = []
    if cds_location is not None:
        constraints.append(
            dc.EnforceTranslation(
                genetic_table=genetic_table, start_codon="keep", location=cds_location
            )
        )
        constraints.append(dc.AvoidChanges(location=(0, cds_location[0])))
        constraints.append(dc.AvoidChanges(location=(cds_location[1], len(seq))))
    # TWIST: What does Twist Bioscience do to improve manufacturability of sequences?
    # TWIST: We won’t introduce any repeat longer than 20 base pairs
    constraints.append(dc.UniquifyAllKmers(20))
    # TWIST: We avoid homopolymer runs of 10 or more bases
    constraints.append(
        dc.AvoidPattern(dc.SequencePattern(r"(.)\1{9,}", is_palyndromic=True))
    )
    # TWIST: We avoid fitted sequences that create global GC% of less than 25% or more than 65% and local GC windows (50 bp) of less than 35% or more than 65%
    constraints.append(dc.EnforceGCContent(0.25, 0.65))
    constraints.append(dc.EnforceGCContent(0.35, 0.65, window=50))
    # TWIST: Once we have optimized your sequence for manufacturing, what does Twist Bioscience do to minimize risk to protein expression of sequences?
    # TWIST: We avoid using rare codons (those with a codon frequency of <8%)
    # TODO: the codon optimization should avoid rare codons as best it can, this hard constraint will occasionally fail, skip for now.
    if cds_location is not None:
        constraints.append(
            dc.AvoidRareCodons(
                0.08, codon_usage_table=aa_to_codons, location=cds_location
            )
        )
    # TWIST: We ensure sequences do not have strong hairpins (ΔG <-8) in the first 48 base pairs of resulting sequences
    # TODO: this is DnaChisel's default config, but is not what Twist uses. skip for now.
    # constraints.append(dc.AvoidHairpins(stem_size=20, hairpin_window=200))
    # TWIST: We check both strands to ensure they don’t contain the enzyme cut sites you asked us to avoid creating
    if avoid_enzymes:
        for enzyme in avoid_enzymes:
            constraints.append(dc.AvoidPattern(f"{enzyme}_site", location=cds_location))
    # TWIST: We avoid introduction of promoter sequences internal to expression sequences by avoiding the creation of strong sigma70 binding sites
    # TODO: this shouldn't be too hard, but we skip for now
    # constraints.append(dc.MotifPssmPattern())
    # TWIST: We avoid sequences that create strong ribosome binding sites (GGAGG and TAAGGAG)
    constraints.append(dc.AvoidPattern(dc.SequencePattern(r"GGAGG|TAAGGAG")))
    # TWIST: We avoid sequences that create terminator sequences (TTTTT or AAAAA)
    constraints.append(dc.AvoidPattern(dc.SequencePattern(r"(A|T)\1{4,}", size=5)))
    # TWIST (error message): >45% of your sequence is composed of small repeats (9bp or longer). This increases complexity. Please break up repeats, perhaps by varying your codon usage.
    # TODO: didn't fix error, so skip for now
    # constraints.append(
    #     dc.AvoidPattern(dc.RepeatedKmerPattern(2, 7))
    # )  # 9-mer repeated twice
    # TODO: don't change codons if not necessary
    if cds_location is not None:
        objectives.append(
            dc.CodonOptimize(
                method=codon_optimization,
                codon_usage_table=aa_to_codons,
                location=cds_location,
            )
        )
    return constraints, objectives


def dnachisel(seq, constraints, objectives):
    if isinstance(seq, Seq):
        seq = str(seq)
    problem = dc.DnaOptimizationProblem(
        sequence=seq, constraints=constraints, objectives=objectives, logger=None
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.to_record(
        with_original_features=True,
        with_sequence_edits=True,
        with_constraints=False,
        with_objectives=False,
    )
