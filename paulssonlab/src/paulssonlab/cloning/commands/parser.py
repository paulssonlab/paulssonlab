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


def unparse_command(cmd, outermost=True):
    # TODO: handle expr args
    if not isinstance(cmd, dict):
        return cmd
    if not outermost and "dest" in cmd and cmd["dest"] is not None:
        dest = f":{cmd['dest']}"
    else:
        dest = ""
    args = ", ".join(
        [unparse_command(arg, outermost=False) for arg in cmd["arguments"]]
    )
    name = cmd["command_name"]
    return f"@{name}{dest}({args})"


# class CommandContext:
#     def __init__(self, registry, context):
#         self.registry = registry
#         self.context = context
#         self.ids = {}
#         self.start_row = {}
#         self.rows = {}

#     def get_id(self, id_, type_):
#         if id_ is None:
#             return None
#         prefix, num = workflow.parse_id(id_)
#         if type_ not in reg.types_for_prefix(prefix):
#             raise ValueError(
#                 f"expecting a collection of type '{type_}' for prefix '{prefix}'"
#             )
#         if num is not None:
#             return workflow.format_id((prefix, num))
#         key = (prefix, type_)
#         id_ = self.ids.get(key)
#         if id_ is None:
#             # TODO
#             # sheet = self.registry.get_sheet((prefix, type_,))
#             # id_, row = workflow.get_next_collection_id(sheet)
#             # id_, row = (("prefix" + prefix, 0), 1)
#             id_, row = self.registry.get_next_id((prefix, type_))
#             self.start_row[key] = row
#         else:
#             id_ = (id_[0], id_[1] + 1)
#         self.ids[key] = id_
#         return workflow.format_id(id_)

#     def insert(self, prefix, type_, row):
#         key = (prefix, type_)
#         if key not in self.rows:
#             self.rows[key] = []
#         self.rows[key].append(row)
#         return row

#     def execute(self):
#         for key, rows in self.rows.items():
#             sheet = registry.get_sheet(key)
#             api.google.sheets.insert_sheet_rows(sheet, self.start_row[key], rows)
#         old_rows = self.rows
#         self.rows = {}
#         return old_rows


# class CommandSemantics(object):
#     def __init__(self, commands, ctx):
#         self.commands = commands
#         self.ctx = ctx

#     def outermost_command(self, ast):
#         if hasattr(ast, "dest") and ast.dest:
#             raise ValueError(
#                 f"cannot specify destination ('{ast.dest}') in outermost command"
#             )
#         else:
#             ast["dest"] = self.ctx.context
#         return self.command(ast)

#     def command(self, ast):
#         if ast.command_name not in self.commands:
#             raise tatsu.semantics.SemanticError(
#                 "command must be one of: {}".format(
#                     ", ".join([f"@{k}" for k in self.commands.keys()])
#                 )
#             )
#         command = self.commands[ast.command_name]
#         if hasattr(ast, "dest"):
#             dest = ast.dest
#         else:
#             dest = None
#         args = ast.arguments
#         res = command(args, dest=dest, ctx=self.ctx)
#         cmd = dict(ast)
#         cmd["arguments"] = [
#             arg["_command"] if isinstance(arg, dict) and "_command" in arg else arg
#             for arg in cmd["arguments"]
#         ]
#         if "_dest" in res and res["_dest"] is not None:
#             cmd["dest"] = res["_dest"]
#         return {"_command": cmd, **res}

#     def int_(self, ast):
#         return int(s)

#     def float_(self, ast):
#         return float(s)

#     def name(self, ast):
#         return {**self.ctx.registry.get(ast.name), "_command": ast.name}
