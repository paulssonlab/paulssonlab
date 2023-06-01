import click

from .inventory.cli import cli as inventory_cli

command_groups = [inventory_cli]

cli = click.CommandCollection(sources=command_groups)
