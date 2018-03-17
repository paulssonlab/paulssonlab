import click
from processing import ingest_nd2_file, quantize_frames
from util import open_zarr_group
import inventory


@click.group()
def cli():
    pass


@cli.command()
@click.argument("in_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("out_path", type=click.Path(writable=True, file_okay=False))
def ingest_nd2(in_path, out_path):
    ingest_nd2_file(in_path, out_path)


@cli.command()
@click.argument("in_path", type=click.Path(exists=True, file_okay=False))
@click.argument("out_path", type=click.Path(writable=True, file_okay=False))
@click.option("--bits", type=int, default=7)
@click.option("--random", type=bool, default=True)
def quantize(in_path, out_path, bits, random):
    in_group = open_zarr_group(in_path)
    out_group = open_zarr_group(out_path)
    quantize_frames(in_group, out_group, bits=bits, random=random)


for command in inventory.commands:
    cli.add_command(command)


def main():
    # logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    # logging.getLogger('urllib3').setLevel(logging.WARNING)
    # logging.getLogger('dropbox').setLevel(logging.WARNING)
    cli()


if __name__ == "__main__":
    main()
