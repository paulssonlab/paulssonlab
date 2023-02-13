import Bio.Restriction
from tatsu.ast import AST
from paulssonlab.cloning.commands.parser import (
    expr_parser,
    expr_list_parser,
    command_parser,
)
from paulssonlab.cloning.sequence import anneal, pcr, assemble
from paulssonlab.cloning.enzyme import re_digest

EXPR_PRIORITY = ["digest", "pcr"]


# def apply_expr(func, *args):
#     # this is a placeholder to wrap all sequence expressions as {"_seq": DsSeqRecord(...)}
#     # until we have a need to pass information beyond the bare sequence
#     args = [a["_seq"] if isinstance(a, dict) and "_seq" in a else a for a in args]
#     return {"_seq": func(*args)}


def eval_exprs_by_priority(s, get_func, priorities=EXPR_PRIORITY):
    if not s.strip():
        return None
    ast = expr_list_parser.parse(s)
    expr = None
    for type_ in priorities:
        priority_expr = [e for e in ast if e["_type"] == type_]
        if len(priority_expr):
            expr = priority_expr[0]
            break
    if expr is None and len(ast):
        expr = ast[0]
    return _eval_command(expr, get_func)


def eval_expr(s, get_func):
    ast = expr_parser.parse(s)
    return _eval_command(ast, get_func)


def _eval_command(ast, get_func):
    if ast is None:
        return None
    type_ = ast["_type"]
    if type_ == "command":
        # TODO: here we can change where we evaluate arguments
        # (i.e., not .get() on the second argument of @Digest, the enzyme name;
        # maybe better just to force quoting the enzyme name)
        ast._frozen = False
        ast.arguments = [_eval_command(a, get_func) for a in ast.arguments]
    if type_ == "name":
        res = get_func(ast["name"])
        return res
    elif type_ == "pcr":
        type_ = "command"
        arguments = [
            _eval_command(ast.template, get_func),
            _eval_command(ast.primer1, get_func),
        ]
        if ast.primer2:
            arguments.append(_eval_command(ast.primer2, get_func))
        ast = AST(command_name="PCR", arguments=arguments)
        # return apply_expr(
        #     pcr,
        #     _eval_expr(expr["template"]),
        #     _eval_expr(expr["primer1"]),
        #     _eval_expr(expr["primer2"]),
        # )
    elif type_ == "digest":
        type_ = "command"
        arguments = [_eval_command(ast.input, get_func), ast.enzyme.name]
        ast = AST(command_name="Digest", arguments=arguments)
        # enzyme_name = expr["enzyme"]["name"]
        # return apply_expr(re_digest_part, _eval_expr(expr["input"]), enzyme)
    elif type_ == "anneal":
        type_ = "command"
        arguments = [
            _eval_command(ast.strand1, get_func),
            _eval_command(ast.strand2, get_func),
        ]
        ast = AST(command_name="Anneal", arguments=arguments)
        # return apply_expr(
        #     anneal, _eval_expr(expr["strand1"]), _eval_expr(expr["strand2"])
        # )
    if type_ == "command":
        command = COMMANDS[ast.command_name]
        # if hasattr(ast, "dest"):
        #     dest = ast.dest
        # else:
        #     dest = None
        args = ast.arguments
        res = command(args)
        # res = command(args, dest=dest, ctx=self.ctx)
        # cmd = dict(ast)
        # TODO: the following is needed for unparsing (after commands pick dest IDs, etc.)
        # cmd["arguments"] = [
        #     arg["_command"] if isinstance(arg, dict) and "_command" in arg else arg
        #     for arg in cmd["arguments"]
        # ]
        # if "_dest" in res and res["_dest"] is not None:
        #     cmd["dest"] = res["_dest"]
        # return {"_command": cmd, **res}
        return res
    else:
        raise ValueError(f"unknown command type '{type_}'")


def eval_command(s, get_func):
    ast = command_parser.parse(s)
    return _eval_command(ast, get_func)


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


def cmd_goldengate(args, dest=None, ctx=None):
    # dest_id = ctx.get_id(dest, "tus")
    seqs = _get_seqs(args)
    product = assemble(seqs, method="goldengate")
    # res = {"_seq": product, "_dest": dest_id}
    res = {"_seq": product}
    return res


def cmd_gib(args, dest=None, ctx=None):
    seqs = _get_seqs(args)
    if len(seqs) == 1:
        # special case to circularize single sequences
        product = seqs[0].assemble(method="gibson")
    else:
        product = assemble(seqs, method="gibson")
    res = {"_seq": product}
    return res


def cmd_pcr(args, dest=None, ctx=None):
    if len(args) == 2:
        primer2 = None
    elif len(args) == 3:
        primer2 = args[2]["_seq"]
    else:
        raise ValueError(
            "@PCR expecting two or three arguments: input, primer1, and optionally primer2"
        )
    template = args[0]["_seq"]
    primer1 = args[1]["_seq"]
    product = pcr(template, primer1, primer2)
    res = {"_seq": product}
    return res


def cmd_digest(args, dest=None, ctx=None):
    if len(args) != 2:
        raise ValueError("@Digest expecting exactly three arguments: input, enzyme")
    input_, enzyme_name = args
    print("!!!", enzyme_name)
    if not hasattr(Bio.Restriction, enzyme_name):
        raise ValueError(f"unknown enzyme '{enzyme_name}'")
    enzyme = getattr(Bio.Restriction, enzyme_name)
    print("input>", input_)
    seqs = _get_seqs(input_)
    products = []
    for seq in seqs:
        products.extend(re_digest(seq, enzyme))
    # product = re_digest_part(input_["_seq"], enzyme)
    product = None  # TODO
    res = {"_seq": product, "_seqs": products}
    return res


def cmd_anneal(args, dest=None, ctx=None):
    if len(args) != 2:
        raise ValueError("@Anneal expecting exactly three arguments: input, enzyme")
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
