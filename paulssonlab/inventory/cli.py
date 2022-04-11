import click


@click.group()
def cli():
    pass


@cli.command()
@click.argument("in_path", type=click.Path(exists=True))
@click.argument(
    "out_path", type=click.Path(writable=True, dir_okay=False), required=True
)
@click.option("--rescan", is_flag=True, default=False)
@click.option(
    "--file-list", type=click.File("r")
)  # TODO: make this a flag, have file list as in_path instead?
@click.option("--skip", multiple=True, default=False)
@click.option("--aggregate", multiple=True, default=["tiff"])
@click.option("--metadata/--no-metadata", default=True)
@click.option("--checksums/--no-checksums", default=True)
@click.option("--duplicates/--no-duplicates", default=False)
def inventory(
    in_path,
    out_path,
    rescan,
    file_list,
    skip,
    aggregate,
    metadata,
    checksums,
    duplicates,
):
    if file_list is not None and rescan:
        raise click.fail("cannot specify both --file-list and --rescan")
    from paulssonlab.inventory.core import make_inventory

    make_inventory(
        in_path,
        out_path,
        rescan=rescan,
        file_list=file_list,
        skip=skip,
        aggregate=aggregate,
        metadata=metadata,
        checksums=checksums,
        duplicates=duplicates,
    )


@cli.command()
@click.argument("file", type=click.Path(exists=True, dir_okay=False))
def inspect_metadata(file):
    from paulssonlab.inventory.core import get_metadata

    # convenience imports
    from paulssonlab.io.metadata import (
        _nikon_tiff_label,
        _nikon_tiff_field,
        parse_nikon_tiff_metadata,
    )

    # import PIL.Image
    import nd2reader
    from IPython import embed

    md = get_metadata(file)
    embed()
