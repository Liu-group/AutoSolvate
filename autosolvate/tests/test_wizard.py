#--------------------------------------------------------------------------------------------------#
# test the input wizard for generating solvated systems
# author: Fangning Ren (2025-12-20)
#--------------------------------------------------------------------------------------------------#
from __future__ import annotations

import os
import pty
import selectors
import signal
import time
from dataclasses import dataclass
from typing import Optional

import shutil
import parmed as pmd
import sys
from pathlib import Path

import pytest

from . import helper_functions as hp


@dataclass
class ReadResult:
    text: str
    timed_out: bool


class PtySession:
    """Minimal PTY session helper for driving interactive CLIs.

    Scope:
    - No environment initialization.
    - No assertions about prompts/correctness.
    - Only supports: start -> read -> send -> close.
    """

    def __init__(self, argv: list[str], *, cwd: Optional[str] = None) -> None:
        self._argv = argv
        self._cwd = cwd
        self._pid: Optional[int] = None
        self._master_fd: Optional[int] = None
        self._buffer = ""

    @property
    def pid(self) -> int:
        if self._pid is None:
            raise RuntimeError("Session not started")
        return self._pid

    def prepare_environment(self) -> None:
        # we should always add the basedir of autosolvate to PYTHONPATH
        # Get the autosolvate package directory
        autosolvate_dir = str(Path(__file__).parent.parent.parent)
        if autosolvate_dir not in sys.path:
            sys.path.insert(0, autosolvate_dir)
        os.environ["PYTHONPATH"] = autosolvate_dir + ":" + os.environ.get("PYTHONPATH", "")

    def start(self) -> None:
        if self._pid is not None:
            raise RuntimeError("Session already started")

        # pty.fork() creates a new child process whose stdin/stdout/stderr are
        # connected to the PTY slave. The parent gets the PTY master fd and can
        # read/write as if it were a terminal.
        pid, master_fd = pty.fork()
        if pid == 0:
            # Child
            self.prepare_environment()
            if self._cwd is not None:
                os.chdir(self._cwd)
            os.execvpe(self._argv[0], self._argv, os.environ.copy())
            raise AssertionError("execvpe failed")

        # Parent
        self._pid = pid
        self._master_fd = master_fd
        # Non-blocking reads keep the event loop responsive. We use selectors to
        # wait until the master fd becomes readable.
        os.set_blocking(master_fd, False)

    def is_alive(self) -> bool:
        if self._pid is None:
            return False
        try:
            pid, status = os.waitpid(self._pid, os.WNOHANG)
        except ChildProcessError:
            return False
        if pid == 0:
            return True
        return False

    def close(self, *, sig: int = signal.SIGTERM, timeout_sec: float = 2.0) -> None:
        if self._pid is None:
            return

        try:
            os.kill(self._pid, sig)
        except ProcessLookupError:
            pass

        deadline = time.time() + timeout_sec
        while time.time() < deadline:
            if not self.is_alive():
                break
            time.sleep(0.05)

        if self.is_alive():
            try:
                os.kill(self._pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

        if self._master_fd is not None:
            try:
                os.close(self._master_fd)
            except OSError:
                pass

        self._pid = None
        self._master_fd = None

    def send(self, line: str) -> None:
        """Send one line (a user turn)."""
        if self._master_fd is None:
            raise RuntimeError("Session not started")
        if not line.endswith("\n"):
            line += "\n"
        os.write(self._master_fd, line.encode("utf-8", errors="replace"))

    def read(self, *, timeout_sec: float = 2.0, max_bytes: int = 64_000) -> ReadResult:
        """Read whatever is available from the child (best-effort)."""
        if self._master_fd is None:
            raise RuntimeError("Session not started")

        sel = selectors.DefaultSelector()
        sel.register(self._master_fd, selectors.EVENT_READ)

        start = time.time()
        read_bytes = 0
        timed_out = False
        try:
            while True:
                if (time.time() - start) >= timeout_sec:
                    timed_out = True
                    break
                if read_bytes >= max_bytes:
                    break
                events = sel.select(timeout=0.2)
                if not events:
                    continue
                try:
                    chunk = os.read(self._master_fd, 4096)
                except BlockingIOError:
                    continue
                except OSError:
                    break
                if not chunk:
                    break
                read_bytes += len(chunk)
                self._buffer += chunk.decode("utf-8", errors="replace")
        finally:
            try:
                sel.unregister(self._master_fd)
            except Exception:
                pass
            sel.close()

        text = self._buffer
        self._buffer = ""
        return ReadResult(text=text, timed_out=timed_out)

class WaitScheduler:
    # a weight time scheduler
    # schedule the wait time based on log change
    # if a log changed, reduce the wait time to 1 second
    # if that does not change, increase the wait time x2 until max wait time 60s is reached
    def __init__(self, initial_wait: float = 1.0, max_wait: float = 60.0) -> None:
        self.current_wait = initial_wait
        self.max_wait = max_wait
        self.last_log_size = 0

    def update(self, log_path: str) -> None:
        try:
            current_size = os.path.getsize(log_path)
        except FileNotFoundError:
            current_size = 0

        if current_size != self.last_log_size:
            self.current_wait = 1.0
            self.last_log_size = current_size
        else:
            self.current_wait = min(self.current_wait * 2, self.max_wait)

def _get_interactive_session() -> PtySession:
    return PtySession(["python", "-m", "autosolvate", "boxgen_interactive"])

def test_wizard_minimal(tmpdir):
    """Minimal smoke test: start wizard, send 'exit', then close."""
    session = _get_interactive_session()
    session.start()
    rr = session.read(timeout_sec=3.0).text 
    assert "welcome" in rr.lower(), rr
    session.send("exit")
    session.read(timeout_sec=3.0)
    session.close()

def test_wizard_boxgen(tmpdir):
    """Test a minimal boxgen wizard session."""
    
    step_inputs = [
        "30",  # Set box size in angstroms
        "1", # choose system type 1 
        "1",  # one solute species
        hp.get_input_dir("naphthalene_neutral.xyz"), # solute file
        "naphthalene",  # solute name
        "1",  # type 1 solute: regular molecule
        "0",  # neutral solute
        "1",  # singlet ground state
        "1",  # single solute solvated in a box
        "no",  # no frcmod or inpcrd file exists
        "yes", # center this solute in the box
        "1",  # single solvent/component
        "2",  # specify solvent amount by density
        "acetonitrile",  # preset solvent name
        "0.776",  # target density in g/cm^3
        "41.05",  # molecular weight in g/mol
        "1.8", # packmol closeness in angstroms
        str(tmpdir),  # output folder
        "yes"   # only write input files, do not run box generation
    ]

    session = _get_interactive_session()
    session.start()
    initial_read = session.read(timeout_sec=2.0)
    # print(initial_read.text)
    for user_input in step_inputs:
        # print(f"Sending input: {user_input}")
        session.send(user_input)
        time.sleep(0.05)  # wait a bit for the process to respond
        read_result = session.read(timeout_sec=0.2)
        # print(read_result.text)
    wait_scheduler = WaitScheduler()
    for i in range(20):
        if os.path.exists("naphthalene-acetonitrile.prmtop") and os.path.exists("naphthalene-acetonitrile.inpcrd") and os.path.exists("naphthalene-acetonitrile.pdb"):
            break
        time.sleep(wait_scheduler.current_wait)
        wait_scheduler.update("autosolvate.log")
    session.close()
    # start to verify generated files
    assert os.path.exists(os.path.join(str(tmpdir), "wizard_input.json")), str(os.listdir(str(tmpdir)))
    assert os.path.exists("naphthalene-acetonitrile.prmtop"), str(os.listdir())
    assert os.path.exists("naphthalene-acetonitrile.inpcrd")
    assert os.path.exists("naphthalene-acetonitrile.pdb")

    prmtop_path = "naphthalene-acetonitrile.prmtop"
    try:
        struct = pmd.load_file(prmtop_path)
    except Exception as e:
        assert False, "Failed to load generated prmtop file." + str(e)

    # Count residues by name
    # struct.residues is a list of Residue objects
    slu_count = sum(1 for res in struct.residues if res.name == 'SLU')
    c3n_count = sum(1 for res in struct.residues if res.name == 'C3N')
    assert slu_count == 1, f"Expected 1 SLU residue, found {slu_count}"
    assert 300 <= c3n_count <= 320, f"Expected 300-320 acetonitrile (C3N), found {c3n_count}"


def test_wizard_perylene_box(tmpdir):
    """Test a pure perylene box generation via the wizard."""
    step_inputs = [
        "55",  # box size in angstroms (Normally 50 Angstrom with 400 molecules matches perylene's crystal density 1.35 g/cm^3, but we use 55 here to reduce the time for packmol to finding the optimal packing configuration)
        "2", # choose system type 2: no solute.
        "1",  # single solvent/component
        "1",  # specify component amount by molecule count
        "custom", # custom component
        hp.get_input_dir("perylene.xyz"),  # to
        "PER",  # residue name
        "no",  # no frcmod or inpcrd file exists
        "400",  # number of solvent molecules
        "2.0", # packmol closeness in angstroms
        str(tmpdir),  # output folder
        "yes"   # only write input files, do not run box generation
    ]

    session = _get_interactive_session()
    session.start()
    initial_read = session.read(timeout_sec=2.0)
    # print(initial_read.text)
    for user_input in step_inputs:
        # print(f"Sending input: {user_input}")
        session.send(user_input)
        time.sleep(0.05)  # wait a bit for the process to respond
        read_result = session.read(timeout_sec=0.2)
        # print(read_result.text)
    wait_scheduler = WaitScheduler()
    for i in range(20):
        if os.path.exists("perylene_box.prmtop") and os.path.exists("perylene_box.inpcrd") and os.path.exists("perylene_box.pdb"):
            break
        time.sleep(wait_scheduler.current_wait)
        wait_scheduler.update("autosolvate.log")
    session.close()
    # start to verify generated files
    assert os.path.exists(os.path.join(str(tmpdir), "wizard_input.json")), str(os.listdir(str(tmpdir)))
    assert os.path.exists("perylene_box.prmtop"), str(os.listdir())
    assert os.path.exists("perylene_box.inpcrd")
    assert os.path.exists("perylene_box.pdb")

    # verify there are 400 perylene molecules in the box
    prmtop_path = "perylene_box.prmtop"
    try:
        struct = pmd.load_file(prmtop_path)
    except Exception as e:
        assert False, "Failed to load generated prmtop file." + str(e)
    perylene_count = sum(1 for res in struct.residues if res.name == 'PER')
    assert perylene_count == 400, f"Expected 400 PER residues, found {perylene_count}"


def _whether_run_ultimate_test() -> bool:
    try:
        hp.get_input_dir("metalcomplex/Fc1N1112_automcpb")
        hp.get_input_dir("IL/TFSI.pdb")
        hp.get_input_dir("EC.xyz")
        hp.get_input_dir("PC.xyz")
        hp.get_input_dir("EMC.xyz")
        return True
    except FileNotFoundError:
        return False

@pytest.mark.skipif(not _whether_run_ultimate_test(), reason="Ultimate complex test skipped due to missing input files.")
def test_wizard_ultimate_task(tmpdir):
    """Test an extremely complex task for the wizard."""
    # This is an ultimate complex task for autosolvate wizard testing.
    # this test requires an:
    # - an trantision metal complex with 2 different legands, only one of which is +1 charged
    # - an counterion with -1 charge and prepared frcmod files from gromacs format
    # - 3 different solvents using weight composition 5:4:1. All of them are custom solvents.
    # ~(づ｡◕‿‿◕｡)づ ENJOY ~

    # TODO: I need to ask my collaborators which orca output files are required for this test. since really running ORCA will take over 3 hours to finish. 

    # before executing this, we have to move a prepared orca result folder to the tmpdir to reduce the time cost.
    orca_result_dir = hp.get_input_dir("metalcomplex/Fc1N1112_automcpb")
    shutil.copytree(orca_result_dir, os.path.join(tmpdir, "Fc1N1112_automcpb"))

    step_inputs = [
        "50", # box size 50 angstrom
        "1",  # choose system type 1
        "2", # there are two different solute species. 

        # solute 1: transition metal complex Fc1N1112
        hp.get_input_dir("metalcomplex/Fc1N1112_automcpb/Fc1N1112.xyz"), # xyz file for solute 1
        "Fc1N1112",  # solute 1 name
        "2", # classify as a transition metal complex
        "+1",  # The Fe ion in this complex is in a +2 oxidation state. One of the ligands is neutral, the other is -1 charged. Overall the complex is +1 charged.
        "1",   # The Fe(II) ion here has all electrons paired.
        "1", # Only one copy of this solute.
        "2", # metal center charge is +2 as Fe(II)
        "skip", # let autosolvate determine the ligand charge by itself.
        "2.8", # Cutoff distance for coordinating ligands
        "yes", # center this solute in the box

        # solute 2: TFSI counterion
        hp.get_input_dir("IL/TFSI.pdb"), # solute 2 structure file
        "TFSI", # solute 2 name
        "1", # regular molecule
        "-1", # -1 charged counterion
        "1", # singlet ground state
        "1", # only one copy of this solute to neutralize the system
        "yes", # the counterion has prepared frcmod and inpcrd files
        "yes", # we have the frcmod file
        hp.get_input_dir("IL/TFSI.frcmod"), # frcmod file path
        hp.get_input_dir("IL/TFSI.mol2"), # provide a mol2 file 

        # MCPB specific inputs
        "orca", # use orca to fit force field parameters
        "/pscratch/sd/f/fren5/orca_6_1_0/orca", # orca executable path. 
        "lanl2dz", # basis set
        "b3lyp", # dft functional
        "16", # single processor
        "4096", # 4096 MB memory each processor
        "yes", # do geometry optimization before fitting
        "skip", # use the detected amberhome.

        # solvent inputs
        "3", # 3 types of solvents
        "2", # use weight portion to specify solvent amount
        "custom", # use custom solvents. 
        hp.get_input_dir("EC.xyz"),
        "EC", # solvent 1 residue name
        "no", # no frcmod or inpcrd file exists
        "1.32", # density in g/cm^3
        "88.06", # molecular weight in g/mol
        "0.4", # weight portion"
        "custom", # solvent 2
        hp.get_input_dir("PC.xyz"),
        "PC", # solvent 2 residue name
        "no", # no frcmod or inpcrd file exists
        "1.2", # density in g/cm^3
        "102.09", # molecular weight in g/mol
        "0.1", # weight portion
        "custom", # solvent 3
        hp.get_input_dir("EMC.xyz"),
        "EMC", # solvent 3 residue name
        "no", # no frcmod or inpcrd file exists
        "1.006", # density in g/cm^3
        "104.15", # molecular weight in g/mol
        "0.5", # weight portion

        # other inputs
        "2.0", # packmol closeness in angstroms
        str(tmpdir),  # output folder
        "yes"   # abort the ultimate complex test to avoid long waiting time.
    ]
    session = _get_interactive_session()
    session.start()
    initial_read = session.read(timeout_sec=2.0)
    # print(initial_read.text)
    for user_input in step_inputs:
        # print(f"Sending input: {user_input}")
        session.send(user_input)
        time.sleep(0.05)  # wait a bit for the process to respond
        read_result = session.read(timeout_sec=0.1)
        # print(read_result.text)
    
    # When tleap generates the topology for this system, there is a bug.
    # Since this system is a ferrocene, one atom has 10 bonds.
    # tleap will report an error when generating the topology file, indicating that the maximum number of bonded partners has been exceeded.
    # Therefore, only check whether the "Fc1N1112-TFSI-EC-PC-EMC.pdb" file is generated.

    waitscheduler = WaitScheduler()
    for i in range(20):
        if os.path.exists("Fc1N1112-TFSI-EC-PC-EMC.pdb"):
            break
        time.sleep(waitscheduler.current_wait)
        waitscheduler.update("autosolvate.log")
    session.close()
    # start to verify generated files
    # print(os.getcwd())
    assert os.path.exists("Fc1N1112-TFSI-EC-PC-EMC.pdb"), str(os.listdir())