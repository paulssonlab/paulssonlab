import tatsu
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

pcr = _type:`pcr` template:dna '<' primer1:name ',' primer2:name '>' ;

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

start = command $ ;

argument
    =
    | quoted_string
    | command
    | float
    | int
    | expr
    ;

command_name = '@' ~ @:/\w+/ ;

command_arglist = '(' ~ ws @+:argument ws {{',' ws @+:argument ws }}* ')' ;

command = _type:`command` command_name:command_name [ ':' dest:/\w+/ ] arguments:command_arglist ;

quoted_string = '"' ~ quoted_string:/[^"]*/ '"' ;

float = float:/\d+\.\d+/ ;

int = int:/\d+/ ;

lookup = '$' ~ name ;

{expr_rules}"""

command_parser = tatsu.compile(command_grammar)


def unparse_expr(expr):
    if isinstance(expr, (list, tuple)):
        return ",".join([unparse_expr(e) for e in expr])
    elif not isinstance(expr, (dict, frozendict)) or "_type" not in expr:
        return expr
    type_ = expr["_type"]
    if type_ == "name":
        return expr["name"]
    elif type_ == "pcr":
        return f"{unparse_expr(expr['template'])}<{unparse_expr(expr['primer1'])},{unparse_expr(expr['primer2'])}>"
    elif type_ == "digest":
        return f"{unparse_expr(expr['input'])}/{unparse_expr(expr['enzyme'])}"
    elif type_ == "anneal":
        return f"{unparse_expr(expr['strand1'])}={unparse_expr(expr['strand2'])}"
    elif type_ == "command":
        if expr.get("dest") is not None:
            dest = f":{expr['dest']}"
        else:
            dest = ""
        args = ", ".join([unparse_expr(arg) for arg in expr["arguments"]])
        name = expr["command_name"]
        return f"@{name}{dest}({args})"
    else:
        raise ValueError(f"unknown expression type: {type_}")
