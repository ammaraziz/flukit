import os
import io
import sys
import shlex
import typer
import shutil
import tempfile
import subprocess

from pandas import DataFrame, read_csv
from rich import print
from pathlib import Path
from signal import Signals
from textwrap import dedent, indent
from importlib_resources import files

try:
    from signal import SIGKILL
except ImportError:
    # A non-POSIX platform
    SIGKILL = None

config_path = files('flukit').joinpath('config/nextclade_datasets')

def update_dataset(lineage: str):
    '''
    nextclade dataset download into config/nextclade_datasets/
    '''
    
    flu_datasets = {
        'h3n2' : 'flu_h3n2_ha',
        'h1n1' : 'flu_h1n1pdm_ha',
        'vic'  : 'flu_vic_ha'
        }
    output_dir = config_path / lineage
    
    # find executable
    nextclade = find_executable(["nextclade"])

    # create call
    call = [
        nextclade, "dataset", "get", 
        "--name", flu_datasets[lineage], 
        "--output-dir", shlex.quote(str(output_dir)),
        "--silent"]

    cmd = " ".join(call)

    # run
    try:
        run_shell_command(cmd, raise_errors = True)
        os.path.isfile(output_dir / "tree.json")
        os.path.isfile(output_dir / "reference.fasta")
        os.path.isfile(output_dir / "sequences.fasta")
        os.path.isfile(output_dir / "genemap.gff")
        print(f"[bold green]Success: Nextclade dataset updated for {lineage}")
    except:
        print("[bold yellow]ERROR: Failed to download Nextclade dataset.[/bold yellow]")
        typer.Exit()

def run_nextclade(
    sequences: Path, 
    lineage: str, 
    output_dir: Path = None
    ) -> DataFrame:
    '''
    run nextclade

    arguments:
        lineage - str       :   h1n1, h3n2, vic
        fasta - Path        :   path to fasta 
        output_dir - Path   :   output directory 
    '''
    nthreads = 2
    dataset = config_path / lineage

    if output_dir is None:
        odir = tempfile.gettempdir()
    else:
        odir = output_dir

    # find executable in path
    nextclade = find_executable(["nextclade"])

    # construct the command syntax
    call = [
        nextclade, "run",
        "--input-dataset", str(dataset), 
        "--output-all", str(odir),
        "--jobs", str(nthreads),
        "--silent", 
        "--in-order",
        str(sequences) 
        ]
    cmd = " ".join(call)

    # run cmd
    try:
        run_shell_command(cmd, raise_errors = True)
        tsv = read_csv(Path(odir) / "nextclade.tsv", sep='\t',  na_filter=False)
        if isinstance(odir, io.BufferedIOBase):
            odir.close()
        return(tsv[['seqName', 'clade', 'qc.overallStatus']])
    except Exception as error:
        print(f"[bold yellow]Calling nextclade failed[bold yellow]: Error: {error}")
        typer.Exit()
    


def run_shell_command(cmd, raise_errors=False, extra_env=None):
    """
    Run the given command string via Bash with error checking.
    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.
    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    #print(f"[bold blue]Running the command:[/bold blue] {cmd}")
    return ShellCommandRunner(cmd, raise_errors=raise_errors, extra_env=extra_env).run()


class ShellCommandRunner:
    def __init__(self, cmd, *, raise_errors=False, extra_env=None):
        self.cmd = cmd
        self.raise_errors = raise_errors
        self.extra_env = extra_env

    def run(self):
        try:
            self.invoke_command()
        except Exception as error:
            self.print_error_message(error)

            if self.raise_errors:
                raise error
            return False

        return True

    def invoke_command(self):
        return subprocess.check_output(
            self.shell_executable + self.shell_args,
            shell=False,
            stderr=subprocess.STDOUT,
            env=self.modified_env,
        )

    @property
    def shell_executable(self):
        if os.name == "posix":
            return ["/bin/bash"]
        else:
            # We try best effort on other systems. For now that means nt/java.
            return ["env", "bash"]

    @property
    def shell_args(self):
        return ["-c", "set -euo pipefail; " + self.cmd]

    @property
    def modified_env(self):
        env = os.environ.copy()

        if self.extra_env:
            env.update(self.extra_env)

        return env

    def print_error_message(self, error):
        if isinstance(error, subprocess.CalledProcessError):
            signal = self.signal_from_error(error)

            if signal:
                message = f"Shell exited from fatal signal {signal.name} when running: {self.cmd}"
            else:
                message = f"Shell exited {error.returncode} when running: {self.cmd}"

            output = (error.output or b'').decode().strip("\n")

            if output:
                message += f"\nCommand output was:\n{indent(output, '  ')}"

            # Bash exits 127 when it cannot find a given command.
            if error.returncode == 127:
                message += "\nAre you sure this program is installed?"

            # Linux's oom-killer issues SIGKILLs to alleviate memory pressure
            elif signal is SIGKILL:
                message += f"\nThe OS may have terminated the command due to an out-of-memory condition."
        elif isinstance(error, FileNotFoundError):
            shell = " and ".join(self.shell_executable)

            message = f"""
                Unable to run shell commands using {shell}!
                Augur requires {shell} to be installed.  Please open an issue on GitHub
                <https://github.com/nextstrain/augur/issues/new> if you need assistance.
                """
        else:
            message = str(error)

        self.print_error(message)

    @staticmethod
    def print_error(message):
        """Prints message to STDERR formatted with textwrap.dedent"""
        print("\nERROR: " + dedent(message).lstrip("\n") + "\n", file=sys.stderr)

    @staticmethod
    def signal_from_error(error):
        """
        Return the :py:class:`signal.Signals` member for the
        :py:attr:`subprocess.CalledProcessError.returncode` of *error*, if any.
        """
        def signal(num):
            try:
                return Signals(num)
            except ValueError:
                return None

        # A grandchild process exited from a signal, which bubbled back up
        # through Bash as 128 + signal number.
        if error.returncode > 128:
            return signal(error.returncode - 128)

        # CalledProcessError documents that fatal signals for the direct child
        # process (bash in our case) are reported as negative exit codes.
        elif error.returncode < 0:
            return signal(abs(error.returncode))

        else:
            return None

def find_executable(names, default = None):
    """
    Return the path to the first executable found in PATH from the given list
    of names.
    Raises a (hopefully helpful) error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of %s in PATH=%s" % (names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception

    return exe