from pathlib import Path
import typer
from snk_cli import CLI
from importlib.metadata import version, PackageNotFoundError

chloroscan = CLI(Path(__file__).parent.parent)

def get_version() -> str:
    try:
        return version("chloroscan")
    except PackageNotFoundError:
        return "unknown"

@chloroscan.app.command()
def version():
    """Print the installed ChloroScan version."""
    print(get_version())

@chloroscan.app.command()
def github():
    """Launch the ChloroScan GitHub page."""
    typer.launch("https://github.com/Andyargueasae/chloroscan")

@chloroscan.app.command()
def docs():
    """Launch the ChloroScan documentation."""
    typer.launch("https://andyargueasae.github.io/chloroscan/")