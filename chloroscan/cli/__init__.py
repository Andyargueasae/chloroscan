from pathlib import Path
import typer
from snk_cli import CLI

chloroscan = CLI(Path(__file__).parent.parent)

@chloroscan.app.command()
def github():
    """ Launch the ChloroScan GitHub page. """
    typer.launch("https://github.com/Andyargueasae/chloroscan")


@chloroscan.app.command()
def docs():
    """ Launch the ChloroScan documentation. """
    typer.launch("https://andyargueasae.github.io/chloroscan/")
