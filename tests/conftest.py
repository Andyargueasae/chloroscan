import re
import difflib
import hashlib
import shutil
import subprocess as sp
from pathlib import Path
from typing import List, NewType, Optional, Union
from appdirs import user_cache_dir
import pytest
import typer
from os import listdir

TargetsType = NewType("TargetsType", Union[str, Path, List[Union[str, Path]]])


class FilesDiffer(Exception):
    pass


class SnakemakePytestException(Exception):
    pass

def get_default_conda_prefix() -> Path:
    return Path(user_cache_dir("chloroscan"))/"conda"


def _md5sum(filename):
    md5 = hashlib.md5()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(128 * md5.block_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def _targets_to_pathlist(pathlist: TargetsType) -> List[Path]:
    if isinstance(pathlist, (Path, str)):
        pathlist = [pathlist]
    for ii, path in enumerate(pathlist):
        if isinstance(path, str):
            pathlist[ii] = Path(path)
    return pathlist


class WorkflowResult:
    def __init__(self, targets: List[Path], work_dir: Path, expected_dir: Path):
        self.targets = targets
        self.work_dir = work_dir
        self.expected_dir = expected_dir

    def get_expected_paths(self, expected_files: Optional[TargetsType] = None) -> List[Path]:
        if expected_files is None:
            return self.targets
        
        return _targets_to_pathlist(expected_files)

    def assert_md5sum(self, checksums:Union[str, List[str]], *, expected_files: Optional[TargetsType] = None,):
        if isinstance(checksums, str):
            checksums = [checksums]

        for expected_file, checksum in zip(self.get_expected_paths(expected_files), checksums):
            generated_path = self.work_dir / expected_file
            generated_checksum = _md5sum(generated_path)
            if generated_checksum != checksum:
                raise SnakemakePytestException(
                    f"The md5 checksum for '{generated_path}' is '{generated_checksum}'."
                    f"This differs from the expected checksum '{checksum}'."
                )

    def assert_contains(self, strings:Union[str, List[str]], *, expected_files: Optional[TargetsType] = None,):
        if isinstance(strings, str):
            strings = [strings]
        
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            text = generated_path.read_text()
            for expected_string in strings:
                if expected_string not in text:
                    raise SnakemakePytestException(
                        f"The file '{generated_path}' does not contain the string '{expected_string}':\n{text}"
                    )       
       
    def assert_not_contains(self, strings:Union[str, List[str]], *, expected_files: Optional[TargetsType] = None,):
        if isinstance(strings, str):
            strings = [strings]
        
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            text = generated_path.read_text()
            for expected_string in strings:
                if expected_string in text:
                    raise SnakemakePytestException(
                        f"The file '{generated_path}' contains the string '{expected_string}':\n{text}"
                    )

    def assert_line_count(self, expected_count=None, min=None, max=None, expected_files: Optional[TargetsType] = None) -> bool:
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            with open(generated_path, "r") as f:
                line_count = sum(1 for _ in f)
                if expected_count:
                    assert expected_count == line_count

                if min:
                    assert line_count >= min

                if max:
                    assert line_count <= max                


    def assert_re(self, patterns:Union[str, List[str]], expected_files: Optional[TargetsType] = None,):
        if isinstance(patterns, str):
            patterns = [patterns]
        
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            text = generated_path.read_text()
            for pattern in patterns:
                if not re.search(pattern, text):
                    raise SnakemakePytestException(
                        f"The file '{generated_path}' does not match with pattern '{pattern}':\n" + 
                        text
                    )

    def assert_exists(self, expected_files: Optional[TargetsType] = None) -> bool:
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            assert generated_path.exists()

    def assert_dir_exists(self, expected_files: Optional[TargetsType] = None) -> bool:
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            assert generated_path.exists()
            assert generated_path.is_dir()

    def assert_glob_count(self, pattern, count=None, min=None, max=None, expected_files: Optional[TargetsType] = None) -> bool:
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            found_count = sum(1 for x in generated_path.glob(pattern))
            if count:
                assert count == found_count

            if min:
                assert found_count >= min

            if max:
                assert found_count <= max                

    def assert_expected(self, expected_files: Optional[TargetsType] = None, diff_on_fail: bool = True) -> bool:
        # Check expected files
        for expected_file in self.get_expected_paths(expected_files):
            generated_path = self.work_dir / expected_file
            expected_path = self.expected_dir / expected_file

            if _md5sum(generated_path) == _md5sum(expected_path):
                continue

            # If different then write the diff
            if diff_on_fail:
                with open(generated_path) as generated, open(expected_path) as expected:

                    generated_lines = generated.readlines()
                    expected_lines = expected.readlines()

                    diff = difflib.unified_diff(generated_lines, expected_lines, n=0)

                    for line in diff:
                        if line.startswith('+'):
                            fg = typer.colors.GREEN
                        elif line.startswith('-'):
                            fg = typer.colors.RED
                        elif line.startswith('^'):
                            fg = typer.colors.BLUE
                        else:
                            fg = None

                        typer.secho(line.rstrip(), fg=fg)

                raise FilesDiffer(f"'{expected_file}' does not match expected.")

        return True


@pytest.fixture
def run_workflow(tmpdir: Path):
    def _run_workflow(targets: TargetsType, *args, expected_dir=None) -> WorkflowResult:
        targets = _targets_to_pathlist(targets)
        work_dir = Path(tmpdir) / "work_dir"
        tests_dir = Path(__file__).parent.resolve()
        expected_dir = expected_dir or tests_dir / "test-data"


        if not expected_dir.exists():
            raise FileNotFoundError(f"Cannot find expected dir '{expected_dir}'.")

        shutil.copytree(
            expected_dir,
            work_dir,\
            ignore=shutil.ignore_patterns('.snakemake'),
            symlinks=True,
        )
        root_dir = Path(__file__).parent.parent
        # snakefile = root_dir/"chloroscan"/"workflow"/"Snakefile"
        # configfile = root_dir/"config"/"ChloroScan.init.yaml"
        input_assembly = root_dir/"tests"/"test-data"/"test_plastid.fasta"
        input_depth_txt = root_dir/"tests"/"test-data"/"contig_depth_test.txt"
        CAT_db = root_dir/"tests"/"test-data"/"CAT_mini"/"db"
        CAT_tax = root_dir/"tests"/"test-data"/"CAT_mini"/"tax"
        input_batch_name = "TEST_OUT"
        output_dir = "TEST_OUT"
        conda_prefix = root_dir/"conda_envs"


        sp.check_output(
            [
                "chloroscan",
                "run",
                "--Inputs-assembly",
                input_assembly,
                "--Inputs-depth-txt",
                input_depth_txt,
                "--Inputs-batch-name",
                input_batch_name,
                "--outputdir",
                output_dir,
                "--cat-database={}".format(CAT_db),
                "--cat-taxonomy={}".format(CAT_tax),
                "-c",
                "32",
                "--use-conda",
                "--force",
                "--reason",
                "-j1",
                "--rerun-triggers=code",
                "--verbose",
                "--latency-wait=15", 
                "--directory={}".format(work_dir),
                "--keep-target-files",
                "--conda-prefix={}".format(conda_prefix),
                *targets,
            ]
        )
        # snakemake -c 11 --latency-wait 15 --use-conda --snakefile ./chloroscan/workflow/Snakefile --configfile config/ChloroScan.init.yaml --directory /home/student.unimelb.edu.au/yuhtong/andy/MMA_organelle_metagenomics --conda-prefix ./conda
        return WorkflowResult(targets, work_dir, expected_dir)

    return _run_workflow