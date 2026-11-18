"""Read ls1-MarDyn *.domainDecompAnalysis files into pandas DataFrames.

The files are written by the DomainDecompAnalyser plugin (src/io/DomainDecompAnalyser).
Changes to the output format there must be reflected in this parser.

Created on 18.11.2026

@author: Georg von Bismarck
"""

import json
import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Final, Optional

import pandas as pd

type Number = int | float
type Point = list[Number]
type Box = list[Point]
type ProcessData = dict[str, list]
type BoxRow = tuple[int, int, Number, Number, Number, Number, Number, Number]
#: Handles one key: value line of a section: (key, raw_value, lineno).
type LineHandler = Callable[[str, str, int], None]

#: Number of constant initial values at the start of every per-simstep list.
N_INITIAL_STEPS: Final[int] = 2

@dataclass(frozen=True)
class DomainDecompGeneralInfo:
    """Contents of the [General information] section.

    Placeholder: the section is currently always empty, so this class has no
    fields yet.
    """


class DomainDecompFormatError(ValueError):
    """Raised when a *.domainDecompAnalysis file does not match the expected format.

    Attributes
    ----------
    path:
        File in which the error was found, if known.
    lineno:
        1-based line number of the error, if it refers to a single line.
    """

    def __init__(
        self, message: str, path: Optional[Path] = None, lineno: Optional[int] = None
    ) -> None:
        self.path: Optional[Path] = path
        self.lineno: Optional[int] = lineno
        location: str = ":".join(str(part) for part in (path, lineno) if part is not None)
        super().__init__(f"{location}: {message}" if location else message)


class DomainDecompWarning(UserWarning):
    """Issued for suspicious but readable data in a *.domainDecompAnalysis file."""


class DomainDecompParser:
    """Parses the lines of one *.domainDecompAnalysis file, section by section.
    """

    #: Per-simstep lists every [Process N] section must contain.
    REQUIRED_STEP_KEYS: Final[tuple[str, ...]] = (
        "traversalTime",
        "averageTraversalTime",
        "numberParticles",
    )
    #: Per-rebuild lists every [Process N] section must contain (may be []).
    REBUILD_KEYS: Final[tuple[str, ...]] = ("rebuildStep", "rebuildDomainBox")
    BOX_COLUMNS: Final[list[str]] = ["xmin", "ymin", "zmin", "xmax", "ymax", "zmax"]
    #: formatVersion values this reader understands.
    SUPPORTED_FORMAT_VERSIONS: Final[frozenset[int]] = frozenset({1})

    def __init__(self, path: Path) -> None:
        self._path: Path = path
        self._has_format_version: bool = False
        self._simsteps: Optional[list[int]] = None
        self._processes: dict[int, ProcessData] = {}
        # Lines before the first section belong to the file header
        self._handle_line: LineHandler = self._parse_header_line

    def parse(self, lines: Iterable[str]) -> tuple[list[int], dict[int, ProcessData]]:
        """Parse all lines and return (simsteps, {rank: {key: list}}). The result is validated
        A Warning is issued if the total numberParticles over all processes is inconsistent.
        """
        for lineno, line in enumerate(lines, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                self._start_section(line[1:-1].strip(), lineno)
                continue
            key, sep, raw = line.partition(":")
            if not sep:
                raise self._error("no ':' found", lineno)
            self._handle_line(key.strip(), raw, lineno)
        return self._validate()

    # --- sections -----------------------------------------------------------

    def _start_section(self, name: str, lineno: int) -> None:
        """Switch the line handler to the section [name]."""
        if name == "General information":
            self._handle_line = self._parse_general_line
        elif name == "Data":
            self._handle_line = self._parse_data_line
        elif name.startswith("Process "):
            self._start_process_section(name, lineno)
        else:
            raise self._error(f"unknown section '[{name}]'", lineno)

    def _start_process_section(self, name: str, lineno: int) -> None:
        """Register the process of a [Process N] header and collect its lines."""
        try:
            rank: int = int(name.removeprefix("Process "))
        except ValueError as exc:
            raise self._error(f"invalid process number in '[{name}]'", lineno) from exc
        if rank in self._processes:
            raise self._error(f"process {rank} appears more than once", lineno)
        data: ProcessData = {}
        self._processes[rank] = data
        self._handle_line = partial(self._parse_process_line, data)

    def _parse_header_line(self, key: str, raw: str, lineno: int) -> None:
        """File header before the first section: only formatVersion."""
        if key != "formatVersion":
            raise self._error(f"'{key}' appears before the first section", lineno)
        try:
            version: int = int(raw.strip())
        except ValueError as exc:
            raise self._error(f"formatVersion '{raw.strip()}' is not an integer", lineno) from exc
        if version not in self.SUPPORTED_FORMAT_VERSIONS:
            raise self._error(
                f"unsupported formatVersion {version}, "
                f"supported: {sorted(self.SUPPORTED_FORMAT_VERSIONS)}",
                lineno,
            )
        self._has_format_version = True

    def _parse_general_line(self, key: str, raw: str, lineno: int) -> None:
        """[General information]: not evaluated yet, see DomainDecompGeneralInfo."""

    def _parse_data_line(self, key: str, raw: str, lineno: int) -> None:
        """[Data]: only simstep is evaluated, unknown keys are skipped with a warning."""
        if key != "simstep":
            warnings.warn(
                f"{self._path}:{lineno}: unexpected key '{key}' in section '[Data]', skipped",
                DomainDecompWarning,
                stacklevel=2,
            )
            return
        self._simsteps = self._parse_list(key, raw, lineno)

    def _parse_process_line(self, data: ProcessData, key: str, raw: str, lineno: int) -> None:
        """[Process N]: every key is stored as a list in data."""
        data[key] = self._parse_list(key, raw, lineno)

    # --- validation ---------------------------------------------------------

    def _validate(self) -> tuple[list[int], dict[int, ProcessData]]:
        """Check that the collected data is complete and consistent and return it."""
        if not self._has_format_version:
            raise self._error("no 'formatVersion' line found")
        if self._simsteps is None:
            raise self._error("no 'simstep' line found")
        if not self._processes:
            raise self._error("no '[Process N]' section found")
        for rank, data in self._processes.items():
            self._validate_required_keys(rank, data)
            self._validate_step_lists(rank, data, len(self._simsteps))
            self._validate_rebuilds(rank, data)
        self._check_particle_count(self._simsteps)
        return self._simsteps, self._processes

    def _validate_required_keys(self, rank: int, data: ProcessData) -> None:
        """Every process must contain all REQUIRED_STEP_KEYS and REBUILD_KEYS."""
        missing: list[str] = [
            key for key in (*self.REQUIRED_STEP_KEYS, *self.REBUILD_KEYS) if key not in data
        ]
        if missing:
            raise self._error(f"process {rank}: missing {', '.join(repr(k) for k in missing)}")

    def _validate_step_lists(self, rank: int, data: ProcessData, n_steps: int) -> None:
        """Every per-simstep list must have one value per simstep."""
        for key, values in data.items():
            if key in self.REBUILD_KEYS:
                continue
            if len(values) != n_steps:
                raise self._error(
                    f"process {rank}: '{key}' has {len(values)} values, "
                    f"expected {n_steps} (length of simstep)"
                )

    def _validate_rebuilds(self, rank: int, data: ProcessData) -> None:
        """One well-formed box per rebuild step."""
        steps: list[int] = data["rebuildStep"]
        boxes: list = data["rebuildDomainBox"]
        if len(steps) != len(boxes):
            raise self._error(
                f"process {rank}: {len(steps)} rebuildStep values "
                f"but {len(boxes)} rebuildDomainBox entries"
            )
        for step, box in zip(steps, boxes):
            if not self._is_box(box):
                raise self._error(
                    f"process {rank}: rebuildDomainBox at rebuildStep {step} is not "
                    f"of the form [[xmin, ymin, zmin], [xmax, ymax, zmax]]"
                )

    def _check_particle_count(self, simsteps: list[int]) -> None:
        """Warn if the total numberParticles over all processes changes."""
        totals: list[int] = [
            sum(counts) for counts in zip(*(data["numberParticles"] for data in self._processes.values()))
        ]
        if not totals:
            return
        expected: int = totals[0]
        deviating: list[int] = [step for step, total in zip(simsteps, totals) if total != expected]
        if deviating:
            warnings.warn(
                f"{self._path}: total numberParticles is not constant: {expected} at simstep "
                f"{simsteps[0]}, but differs in {len(deviating)} simstep(s), first at simstep "
                f"{deviating[0]} (min {min(totals)}, max {max(totals)})",
                DomainDecompWarning,
                stacklevel=2,
            )

    @staticmethod
    def _is_box(box: object) -> bool:
        """True if box is [[x, y, z], [x, y, z]]."""
        return (
            isinstance(box, list)
            and len(box) == 2
            and all(isinstance(point, list) and len(point) == 3 for point in box)
        )

    # --- helpers ------------------------------------------------------------

    def _parse_list(self, key: str, raw: str, lineno: int) -> list:
        """Parse the value after key: as a JSON list."""
        try:
            value: object = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise self._error(f"value of '{key}' is not a valid list", lineno) from exc
        if not isinstance(value, list):
            raise self._error(f"value of '{key}' is not a list", lineno)
        return value

    def _error(self, message: str, lineno: Optional[int] = None) -> DomainDecompFormatError:
        """Create a format error for this file."""
        return DomainDecompFormatError(message, self._path, lineno)


def _build_steps(simsteps: list[int], processes: dict[int, ProcessData], skip_initial: bool) -> pd.DataFrame:
    """Measurements per simstep; expects data validated by DomainDecompParser."""
    start: int = N_INITIAL_STEPS if skip_initial else 0
    frames: list[pd.DataFrame] = []
    for rank, data in processes.items():
        metrics: dict[str, list[Number]] = {
            k: v for k, v in data.items() if k not in DomainDecompParser.REBUILD_KEYS
        }
        index: pd.MultiIndex = pd.MultiIndex.from_product(
            [[rank], simsteps[start:]], names=["process", "simstep"]
        )
        frames.append(pd.DataFrame({k: v[start:] for k, v in metrics.items()}, index=index))
    return pd.concat(frames).sort_index()


def _build_domains(processes: dict[int, ProcessData]) -> Optional[pd.DataFrame]:
    """Domain boxes per rebuild; expects data validated by DomainDecompParser."""
    rows: list[BoxRow] = []
    for rank, data in processes.items():
        steps: list[int] = data["rebuildStep"]
        boxes: list[Box] = data["rebuildDomainBox"]
        for step, ((xmin, ymin, zmin), (xmax, ymax, zmax)) in zip(steps, boxes, strict=True):
            rows.append((rank, step, xmin, ymin, zmin, xmax, ymax, zmax))
    if not rows:
        return None
    return (
        pd.DataFrame(rows, columns=["process", "rebuildStep", *DomainDecompParser.BOX_COLUMNS])
        .set_index(["process", "rebuildStep"])
        .sort_index()
    )


def read_domain_decomp(path: str | Path, skip_initial: bool = False) -> tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[DomainDecompGeneralInfo]]:
    """Read a *.domainDecompAnalysis file.

    Parameters
    ----------
    path:
        Path to the file.
    skip_initial:
        If True, the first N_INITIAL_STEPS (=2) simsteps are
        dropped, because their values are constant initial values.

    Returns
    -------
    steps:
        Measurements per simstep, MultiIndex (process, simstep); columns e.g.
        traversalTime, averageTraversalTime, numberParticles.
    domains:
        Domain boxes per rebuild, MultiIndex (process, rebuildStep);
        columns xmin, ymin, zmin, xmax, ymax, zmax. None if
        rebuildStep is empty for all processes.
    general_info:
        Contents of the [General information] section. Currently always
        None, because the section holds no information yet.

    Raises
    ------
    DomainDecompFormatError
        If the file does not match the expected format.
    """
    file_path: Path = Path(path)
    with file_path.open(encoding="utf-8") as f:
        simsteps, processes = DomainDecompParser(file_path).parse(f)
    general_info: Optional[DomainDecompGeneralInfo] = None
    return (
        _build_steps(simsteps, processes, skip_initial),
        _build_domains(processes),
        general_info,
    )


if __name__ == "__main__":
    import sys

    steps, domains, general_info = read_domain_decomp(
        sys.argv[1], skip_initial="--skip-initial" in sys.argv
    )
    print(general_info)
    print(steps.info())
    print(steps.head())
    print(domains if domains is None else domains.head())
