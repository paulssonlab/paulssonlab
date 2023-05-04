import tatsu
from tatsu.ast import AST
from frozendict import frozendict

expr_rules = r"""ws = /\s*/ ;

expr_list = ws @+:expr ws {',' ws @+:expr ws }* ;

expr
    =
    | restriction_digest
    | pcr
    | anneal
    | name
    ;

name = _type:`name` name:/[A-Za-z0-9\._-]+/ ;

anneal = _type:`anneal` strand1:name '=' strand2:name ;

dna = anneal | name ;

pcr = _type:`pcr` template:dna '<' primer1:name [ ',' primer2:name ] '>' ;

restriction_digest = _type:`digest` input:expr '/' enzyme:name ;"""

expr_grammar = rf"""@@grammar::Expr
@@whitespace :: //

start = expr $ ;

{expr_rules}
"""

expr_parser = tatsu.compile(expr_grammar)

expr_list_grammar = rf"""@@grammar::Expr
@@whitespace :: //

start = expr_list $ ;

{expr_rules}
"""

expr_list_parser = tatsu.compile(expr_list_grammar)

command_grammar = rf"""@@grammar::Command
@@whitespace :: //

start = command | expr $ ;

arg
    =
    | quoted_string
    | command
    | float
    | int
    | expr
    ;

command_name = '@' ~ @:/\w+/ ;

command_arglist = '(' ~ ws @+:arg ws {{',' ws @+:arg ws }}* ')' ;

command = _type:`command` command:command_name args:command_arglist ;

quoted_string = _type:`quoted_string` '"' ~ quoted_string:/[^"]*/ '"' ;

float = _type:`float` float:/\d+\.\d+/ ;

int = _type:`int` int:/\d+/ ;

{expr_rules}"""

command_parser = tatsu.compile(command_grammar)


class Name(str):
    def __repr__(self):
        return f"name:'{self}'"


def normalize_ast(ast, recursive=True):
    if recursive:
        _normalize_ast = normalize_ast
    else:
        _normalize_ast = lambda x: x
    if not isinstance(ast, (dict, AST)):
        return ast
    type_ = ast["_type"]
    if type_ == "command":
        return dict(
            _type="command",
            command=ast["command"],
            args=[_normalize_ast(a) for a in ast["args"]],
            kwargs=ast.get("kwargs") or {},
        )
    elif type_ == "name":
        return Name(ast["name"])
    elif type_ == "quoted_string":
        return ast["quoted_string"]
    elif type_ == "int":
        return int(ast["int"])
    elif type_ == "float":
        return float(ast["float"])
    elif type_ == "pcr":
        type_ = "command"
        args = [
            _normalize_ast(ast.template),
            _normalize_ast(ast.primer1),
        ]
        if ast.primer2:
            args.append(_normalize_ast(ast.primer2))
        return dict(_type="command", command="PCR", args=args)
    elif type_ == "digest":
        type_ = "command"
        args = [_normalize_ast(ast.input), ast.enzyme.name]
        return dict(_type="command", command="Digest", args=args)
    elif type_ == "anneal":
        type_ = "command"
        args = [
            _normalize_ast(ast.strand1),
            _normalize_ast(ast.strand2),
        ]
        return dict(_type="command", command="Anneal", args=args)
    else:
        raise ValueError(f"unknown ast type: {type_}")


def unparse_expr(expr):
    if isinstance(expr, (list, tuple)):
        return ",".join([unparse_expr(e) for e in expr])
    elif not isinstance(expr, (dict, frozendict)) or "_type" not in expr:
        return expr
    type_ = expr["_type"]
    if type_ == "name":
        return expr["name"]
    elif type_ == "quoted_string":
        return f"\"{expr['quoted_string']}\""
    elif type_ == "int":
        return expr["int"]
    elif type_ == "float":
        return expr["float"]
    elif type_ == "pcr":
        primers = unparse_expr(expr["primer1"])
        if expr.get("primer2"):
            primers += f",{unparse_expr(expr['primer2'])}"
        return f"{unparse_expr(expr['template'])}<{primers}>"
    elif type_ == "digest":
        return f"{unparse_expr(expr['input'])}/{unparse_expr(expr['enzyme'])}"
    elif type_ == "anneal":
        return f"{unparse_expr(expr['strand1'])}={unparse_expr(expr['strand2'])}"
    elif type_ == "command":
        args = ", ".join([unparse_expr(arg) for arg in expr["args"]])
        name = expr["command"]
        return f"@{name}({args})"
    else:
        raise ValueError(f"unknown expression type: {type_}")
