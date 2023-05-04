import Bio.Restriction
from tatsu.ast import AST
from itertools import repeat
from paulssonlab.cloning.commands.parser import (
    expr_parser,
    expr_list_parser,
    command_parser,
    normalize_ast,
    Name,
)
from paulssonlab.cloning.sequence import anneal, pcr, assemble
from paulssonlab.cloning.enzyme import re_digest
from paulssonlab.cloning.workflow import re_digest_part


def eval_ast(ast, ctx=None):
    ast = normalize_ast(ast, recursive=False)
    if ast.__class__ == Name:
        return ctx["registry"].get(str(ast))
    type_ = ast.get("_type")
    if type_ == "command":
        command = COMMANDS[ast["command"]]
        args = ast.get("args") or []
        kwargs = ast.get("kwargs") or {}
        res = command(args, kwargs, ctx=ctx)
        return res
    else:
        raise ValueError(f"unknown command type '{type_}'")


def _eval_args(args, eval_args, kwargs={}, eval_kwargs=set(), ctx=None):
    args = [
        eval_ast(arg, ctx=ctx) if eval_ else arg for eval_, arg in zip(eval_args, args)
    ]
    kwargs = {
        k: eval_ast(v, ctx=ctx) if k in eval_kwargs else v for k, v in kwargs.items()
    }
    return args, kwargs


# TODO
# def eval_expr(s, get_func):
#     ast = expr_parser.parse(s)
#     return eval_ast(ast, get_func)

# TODO
# def eval_command(s, get_func):
#     ast = command_parser.parse(s)
#     return eval_ast(ast, get_func)


def _get_seqs(args):
    if not isinstance(args, list):
        args = [args]
    seqs = []
    for arg in args:
        if "_seqs" in arg:
            seqs.extend(arg["_seqs"])
        else:
            seqs.append(arg["_seq"])
    return seqs


def cmd_goldengate(args, kwargs, ctx=None):
    args, _ = _eval_args(args, repeat(True), ctx=ctx)
    seqs = _get_seqs(args)
    product = assemble(seqs, method="goldengate")
    res = {"_seq": product}
    return res


def cmd_gib(args, kwargs, ctx=None):
    args, _ = _eval_args(args, repeat(True), ctx=ctx)
    seqs = _get_seqs(args)
    if len(seqs) == 1:
        # special case to circularize single sequences
        product = seqs[0].assemble(method="gibson")
    else:
        product = assemble(seqs, method="gibson")
    res = {"_seq": product}
    return res


def cmd_pcr(args, kwargs, ctx=None):
    args, _ = _eval_args(args, repeat(True), ctx=ctx)
    if len(args) == 2:
        primer2 = None
    elif len(args) == 3:
        primer2 = args[2]["_seq"]
    else:
        raise ValueError(
            "@PCR expecting two or three args: input, primer1, and optionally primer2"
        )
    template = args[0]["_seq"]
    primer1 = args[1]["_seq"]
    product = pcr(template, primer1, primer2)
    res = {"_seq": product}
    return res


def cmd_digest(args, kwargs, ctx=None):
    if len(args) != 2:
        raise ValueError("@Digest expecting exactly two args: input, enzyme")
    args, _ = _eval_args(args, [True, False], ctx=ctx)
    input_, enzyme_name = args
    if not hasattr(Bio.Restriction, enzyme_name):
        raise ValueError(f"unknown enzyme '{enzyme_name}'")
    enzyme = getattr(Bio.Restriction, enzyme_name)
    seqs = _get_seqs(input_)
    products = []
    # use re_digest_part if we have overhangs in kwargs (raise exception if have more than one seq)
    if kwargs.get("overhangs"):
        if len(seqs) != 1:
            raise ValueError("cannot digest multiple sequences if overhangs are given")
        product = re_digest_part(seqs[0], enzyme, overhangs=kwargs["overhangs"])
    else:
        product = None
        for seq in seqs:
            products.extend(re_digest(seq, enzyme))
    res = {"_seq": product, "_seqs": products}
    return res


def cmd_anneal(args, kwargs, ctx=None):
    if len(args) != 2:
        raise ValueError("@Anneal expecting exactly two args: strand1, strand2")
    args, _ = _eval_args(args, repeat(True), ctx=ctx)
    strand1, strand2 = args
    product = anneal(strand1["_seq"], strand2["_seq"])
    res = {"_seq": product}
    return res


COMMANDS = {
    "GG": cmd_goldengate,
    "Gib": cmd_gib,
    "PCR": cmd_pcr,
    "Digest": cmd_digest,
    "Anneal": cmd_anneal,
}
