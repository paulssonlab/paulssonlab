import tatsu

expr_rules = r"""ws = /\s*/ ;

expr_list = ws @+:expr ws {',' ws @+:expr ws }* ;

expr
    =
    | restriction_digest
    | pcr
    | anneal
    | name
    ;

name = _type:`name` name:/[A-Za-z0-9\.]+/ ;

anneal = _type:`anneal` strand1:name '=' strand2:name ;

dna = anneal | name ;

pcr = _type:`pcr` template:dna '<' primer1:name ',' primer2:name '>' ;

restriction_digest = _type:`digest` input:expr '/' enzyme:name ;"""

expr_grammar = rf"""@@grammar::Expr
@@whitespace :: //

start = expr_list $ ;

{expr_rules}
"""

expr_parser = tatsu.compile(expr_grammar)

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

command = command_name:command_name arguments:command_arglist ;

quoted_string = '"' ~ quoted_string:/[^"]*/ '"' ;

float = float:/\d+\.\d+/ ;

int = int:/\d+/ ;

lookup = '$' ~ name ;

{expr_rules}"""

command_parser = tatsu.compile(command_grammar)
