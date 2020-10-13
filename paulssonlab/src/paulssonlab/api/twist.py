import requests
import time
from paulssonlab.cloning.sequence import get_seq
from paulssonlab.cloning.util import enzymes_to_names

API_URL = "https://twist-api.twistbioscience-staging.com"

OPTIMIZATION_CONCURRENCY = 30

SCORING_ERRORS = """4000 - General error
4001 - Problematic sequence
4002 - Repeats or extreme high/low GC in the highlighted region may have caused this problem
4003 - Provided sequence is invalid
4004 - Provided sequence contains invalid character(s)
4100 - Invalid sequence length
4101 - Sequence is too short
4102 - Sequence is too long
4103 - Sequences longer than 1,700 bp elevate risk marginally.  Consider splitting your sequence into smaller pieces to increase the likelihood of success
4104 - Sequences longer than 3,100 bp elevate risk.  Consider splitting your sequence into smaller pieces to increase the likelihood of success
4200 - Invalid GC content
4201 - The overall GC content of your sequence must be under 65%. Under 60% will be optimal for success
4202 - The overall GC content of your sequence must be over 25%
4203 - The difference in GC content between the highest-GC and lowest-GC 50bp windows exceeds 52%. Please even out the peaks/troughs in GC content of your sequence. Lowering the peaks will be optimal for success
4300 - Secondary structure
4301 - Hairpin detected
4303 - Long direct repeat (greater or equal to 20bp) detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success
4304 - Direct Repeat with high Tm (greather than 60C) detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success
4305 - More than 45% of your sequence is composed of small repeats (9bp or longer). This increases complexity. Please break up repeats, perhaps by varying your codon usage
4306 - Long, low homology repeat region detected. Please break up repeats. Fewer/lower homology repeats will be optimal for success.
4401 - Attempted hierarchical design but could not find an acceptable split point
4402 - We are unable to make the sequence as is. Please try to run our Codon Optimization or try to split the gene into two parts
4403 - Attempted hierarchical design but could not find an acceptable split point
4404 - Design request results in fragment sizes less than 300bp; please retry design with fewer fragments
4405 - Design request results in fragment sizes greater than 1,800bp; please retry design with more fragments
4501 - His tags consisting of 5 or more identical codons increase complexity. Please vary the codons used (e.g. CAT CAC CAT... instead of CAT CAT CAT...)
4502 - CpG multimeric segments of 14 or more bases increase complexity. Please break up these low-complexity sequences
4503 - Long homopolymer stretches increase complexity. Please break up long homopolymers
4504 - Sequence contains one or more Gateway cloning att sites
4505 - Sequence contains one or more impermissible sequences that complicate manufacture or QC.
4506 - Clonal gene contains sub-sequence with high homology to CCDB.
4507 - Sequence contains one or more sequences that may pose specific risk to downstream experiments.
4508 - Unable to design primers for this construct.  Please increase the GC content in the first and last 60 bases of your sequence.
5001 - Hierarchical design exception
5002 - Fragment design exception"""
SCORING_ERRORS = dict([line.split(" - ") for line in SCORING_ERRORS.split("\n")])


def _format_seq(seq):
    return str(get_seq(seq)).upper()


def auth_headers(config):
    return {
        "Authorization": "JWT {}".format(config["api_token"]),
        "X-End-User-Token": config["end_user_token"],
    }


def score(config, seqs, gene_type="NON_CLONED_GENE"):
    data = [
        {
            "sequences": [_format_seq(seq)],
            "name": name,
            "external_id": name,
            "type": gene_type,
        }
        for name, seq in seqs.items()
    ]
    scoring_jobs = _score(config, data)
    # return constructs
    ids = [j["id"] for j in scoring_jobs]
    scored_constructs = _poll_score(config, ids)
    scored_constructs = {s["external_id"]: s for s in scored_constructs}
    # keep same order as input
    scored_constructs = {name: scored_constructs[name] for name in seqs.keys()}
    return scored_constructs


def _score(config, data):
    r = requests.post(
        f"{API_URL}/v1/users/{config['ecommerce_username']}/constructs/",
        json=data,
        headers=auth_headers(config),
    )
    r.raise_for_status()
    return r.json()


# TODO: async version
def _poll_score(config, ids, delay=5, max_retries=20):
    # TODO: this probably maxes out at 500 returned constructs (inferred from swagger docs example)
    # TODO: this is also severely limited by GET URL param length limits, need to ask Twist to add a POST method
    for i in range(0, max_retries):
        r = requests.get(
            f"{API_URL}/v1/users/{config['ecommerce_username']}/constructs/describe/",
            params={"id__in": ",".join(ids), "scored": True},
            headers=auth_headers(config),
        )
        if r.status_code == 200:
            scored_constructs = r.json()
            if set(c["id"] for c in scored_constructs if c["scored"] is True) == set(
                ids
            ):
                return scored_constructs
        time.sleep(delay)
    raise Exception("Twist scoring timed out")


def optimize_codons(config, seqs, avoid_enzymes=None, organism="Escherichia coli"):
    # TODO: endpoint can only handle max batches of 30
    data = [
        {
            "sequence": _format_seq(seq[0]),
            "external_id": name,
            "optimization_start": seq[1][0],
            "optimization_len": seq[1][1] - seq[1][0],
            "avoid_introducing": enzymes_to_names(avoid_enzymes)
            if avoid_enzymes
            else [],
            "organism": organism,
        }
        for name, seq in seqs.items()
    ]
    optimization_jobs = _optimize_codons(config, data)
    ids = [j["id"] for j in optimization_jobs]
    optimized_constructs = _poll_optimization(config, ids)
    optimized_constructs = {o["external_id"]: o for o in optimized_constructs}
    # keep same order as input
    optimized_constructs = {name: optimized_constructs[name] for name in seqs.keys()}
    return optimized_constructs


def _optimize_codons(config, data):
    r = requests.post(
        f"{API_URL}/v1/users/{config['ecommerce_username']}/codon-optimizations/",
        json=data,
        headers=auth_headers(config),
    )
    r.raise_for_status()
    return r.json()


# TODO: async version
def _poll_optimization(config, ids, delay=5, max_retries=20):
    # TODO: this probably maxes out at 500 returned constructs (inferred from swagger docs example)
    # TODO: this is also severely limited by GET URL param length limits, need to ask Twist to add a POST method
    for i in range(0, max_retries):
        r = requests.get(
            f"{API_URL}/v1/users/{config['ecommerce_username']}/codon-optimizations/",
            params={"id__in": ",".join(ids), "completed": True},
            headers=auth_headers(config),
        )
        if r.status_code == 200:
            optimization_jobs = r.json()
            completed_ids = set(
                j["id"] for j in optimization_jobs if j["completed"] is True
            )
            if completed_ids == set(ids):
                return optimization_jobs
        time.sleep(delay)
    raise Exception("Twist codon optimization timed out")


def score_and_optimize(config, seqs, avoid_enzymes=None, organism="Escherichia coli"):
    scored_constructs = score(config, {name: seq[0] for name, seq in seqs.items()})
    need_to_optimize = {
        name: seq
        for name, seq in seqs.items()
        if scored_constructs[name]["score"]
        != "STANDARD"  # TODO: need production API to test this
    }
    optimized_constructs = optimize_codons(
        config, need_to_optimize, avoid_enzymes=avoid_enzymes, organism=organism
    )
    # TODO: simplify output?
    return optimized_constructs
