"""Common model classes used across all workflows."""
from dataclasses import asdict, dataclass, field
from enum import Enum
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..util import get_named_logger  # noqa: ABS101

logger = get_named_logger("Models")


@dataclass
class WorkflowBaseModel:
    """Common things for stuff in the model."""

    def get(
        self,
        field_name: str,
        title: bool = True,
        **kwargs
    ):
        """Get reportable field tuple."""
        field_info = self.__dataclass_fields__.get(field_name)
        # provide an empty string default title to minimise drama
        field_title = field_info.metadata.get("title", "")
        value = self.get_reportable_value(field_name=field_name, **kwargs)
        if title:
            return (field_title, value)
        return value

    def get_reportable_value(
            self,
            field_name: str,
            *,
            decimal_places: int = None,
            default_value: str = "N/A") -> Optional[str]:
        """Get the value of a value and make it reportable."""
        # Get the field info using the field name
        field_info = self.__dataclass_fields__.get(field_name)
        if field_info is None:
            raise AttributeError(
                f"{field_name!r} is not a field on {self.__class__.__name__}"
            )

        value = getattr(self, field_name)

        if value is None:
            return default_value

        if isinstance(value, (int, float)):
            if decimal_places:
                value = round(value, decimal_places)
            if value < 0.0001 or value > 99999999:
                value = f"{value:.2E}"
        else:
            if decimal_places:
                raise TypeError(
                    "decimal_places is not a supported argument for a non-numeric.")

        unit = field_info.metadata.get('unit')

        if unit:
            return f"{value} {unit}"

        return str(value)


class SampleType(str, Enum):
    """The type of the sample."""

    no_template_control = "no_template_control"
    positive_control = "positive_control"
    negative_control = "negative_control"
    test_sample = "test_sample"

    def friendly_name(self):
        """Convert sample type to string."""
        return self.name.replace("_", " ").capitalize()


@dataclass
class SampleIdentifier:
    """Additional identifiers for a sample."""

    name: str = field(
        metadata={
            "title": "Identifier name",
            "Description": "The name of the sample identifier"})
    value: str = field(
        metadata={
            "title": "Identifier value",
            "Description": "The value of the sample identifier"})


@dataclass
class CheckResult:
    """
    A result of some check the workflow has performed.

    This can be at sample or workflow level.
    """

    check_category: str = field(
        metadata={
            "title": "Check category",
            "description": "The category of the check"})
    check_name: str = field(
        metadata={
            "title": "Check name",
            "description": "The name of the check"})
    check_pass: bool = field(
        metadata={
            "title": "Check pass",
            "description": "If true the check has passed"})
    check_threshold: str | None = field(
        default=None, metadata={
            "title": "Check threshold",
            "description": "The threshold for the check, useful for reporting later"})

    categories = {}

    def friendly_check_category(self):
        """Convert category to string."""
        if self.check_category not in self.categories:
            raise ValueError(f"{self.check_category} has no friendly name")
        return self.categories[self.check_category]

    def friendly_check_name(self):
        """Convert check name to string."""
        return self.check_name.replace("_", " ").capitalize()


@dataclass
class ResultsContents:
    """Placeholder class for results contents."""

    pass


@dataclass
class Sample:
    """A sample sheet entry and its corresponding checks and related results."""

    alias: str = field(
        metadata={
            "title": "Sample alias",
            "description": "The alias for the sample given by the user"})
    sample_type: SampleType = field(
        metadata={
            "title": "Sample type",
            "description": "The type of the sample"})
    sample_pass: bool = field(
        metadata={
            "title": "Sample pass",
            "description": "If true the sample has passed workflow checks"})
    barcode: str | None = field(
        default=None,
        metadata={
            "title": "Sample barcode",
            "description": "The physical barcode assigned to the sample"})
    additional_identifiers: List[SampleIdentifier] = field(
        default_factory=list, metadata={
            "title": "Additional sample identifiers",
            "description": "Additional identifiers for the sample"})
    sample_checks: list[CheckResult] = field(
        default_factory=list, metadata={
            "title": "Sample checks",
            "description": "An array of checks performed on the sample"})
    results: ResultsContents | None = field(
        default=None, metadata={
            "title": "Sample results",
            "description": "Further specific workflow results for this sample"})
    config:  Dict[str, Any] | None = field(
        default=None, metadata={
            "title": "Sample configuration",
            "description": """Sample specific config parameters
            used for running analysis"""})

    def __post_init__(self):
        """Determine overall status for a sample given the individual check results."""
        self.sample_pass = all(
            check.check_pass for check in self.sample_checks)

    def get_sample_identifier(self, sample_identifier):
        """Get a sample identifier given the identifier name."""
        for identifier in self.additional_identifiers:
            if identifier.name == sample_identifier:
                return identifier.value
        raise KeyError("Sample identifier not found")

    def set_sample_identifier(self, name, value):
        """Set a sample identifier."""
        sample_identifier = SampleIdentifier(
            name=name,
            value=value)
        self.additional_identifiers.append(sample_identifier)
        return self.additional_identifiers

    def to_json(self, filename):
        """Save class as JSON."""
        with open(filename, 'w') as f:
            json.dump(asdict(self), f, default=str, indent=2)

    def get_reportable_qc_status(self, max_criteria=4):
        """Store global status of the sample and list of QC criteria to show.

        :params max_criteria: Maximum number of criteria to be reported.
        """
        # Store global status: pass/ failed
        qc_global_status = {"status": self.sample_pass, "scope": "QC status"}
        qc_criteria = []
        if self.sample_pass:
            qc_criteria.append(
                    {"status": self.sample_pass, "scope": "All acceptance criteria met"}
            )
        else:
            # Report failed criteria until a maximum value
            for qc in self.sample_checks:
                if not qc.check_pass:  # append criteria if failed
                    qc_criteria.append(
                        {
                            "status": qc.check_pass,
                            "category": qc.friendly_check_category(),
                            "scope": qc.friendly_check_name(),
                        }
                    )
            if len(qc_criteria) > max_criteria:
                # Replace all the failed criteria, with a sentence with the number
                # instead of listing all of them.
                # Set status to False as more than max_criteria are failed.
                qc_criteria = [
                    {
                        "status": False,
                        "scope": f"{len(qc_criteria)} acceptance criteria",
                    },
                ]
        return qc_global_status, qc_criteria


@dataclass
class RunStats:
    """Basic run statistics for the entire run."""

    total_reads: int | None = field(
        default=None, metadata={
            "title": "Total reads",
            "description": "Total number of reads on run"})
    total_ambiguous_reads: int | None = field(
        default=None, metadata={
            "title": "Total ambiguous reads",
            "description": "Number of reads of unknown provenance"})
    total_unaligned_reads: int | None = field(
        default=None, metadata={
            "title": "Total unaligned reads",
            "description": "Number of unaligned reads"})


@dataclass
class WorkflowResult(WorkflowBaseModel):
    """
    Definition for results that will be returned by this workflow.

    This structure will be passed through by Gizmo speaking clients
    as WorkflowInstance.results.
    """

    samples: list[Sample] = field(
        metadata={
            "title": "Samples",
            "description": "Samples in this workflow instance"})
    workflow_pass: bool | None = field(
        default=None, metadata={
            "title": "Workflow pass",
            "description": "True if this workflow instance passes all checks"})
    workflow_checks: list[CheckResult] = field(
        default_factory=list, metadata={
            "title": "Workflow checks",
            "description": "An array of checks performed on the workflow instance"})
    run_stats: RunStats | None = field(
        default=None, metadata={
            "title": "Samples",
            "description": "Basic run statistics"})
    client_fields: dict[str, Any] | None = field(
        default_factory=dict, metadata={
            "title": "Client fields",
            "description": "Arbitrary key-value pairs provided by the client"})
    versions: dict[str, Any] | None = field(
        default_factory=dict, metadata={
            "title": "Analysis tool versions",
            "description": """Key-value pairs collecting the
            software used and the corresponding versions"""})
    params: dict[str, Any] | None = field(
        default_factory=dict, metadata={
            "title": "Pertinent parameters",
            "description": """Key-value pairs with the
            options chosen by the user"""})

    def load_client_fields(self, filename):
        """Load client fields."""
        with open(filename) as f:
            try:
                client_fields = json.loads(f.read())
                # convert any lists into strings for display
                for key, value in client_fields.items():
                    if isinstance(value, list):
                        client_fields[key] = ', '.join(value)
            except json.decoder.JSONDecodeError:
                client_fields = {"error": "Error parsing client fields file."}

        self.client_fields = client_fields
        return self.client_fields

    def load_params(self, params_json, keep=None):
        """Create a workflow params dict."""
        params_json = Path(params_json)
        if keep is None:
            keep = []
        if not params_json.is_file():
            raise FileNotFoundError(f"No such file: {params_json}")
        with open(params_json, "r") as f:
            try:
                params_dict = json.loads(f.read())
                self.params = {
                    k: v for k, v in params_dict.items() if k in set(keep)
                }
                return self.params
            except ValueError:
                raise ValueError(f"Invalid JSON file: {params_json}")

    def load_versions(self, versions_path):
        """Create a version list of dict."""
        versions_path = Path(versions_path)
        if not versions_path.exists():
            raise FileNotFoundError(f"No such file: {versions_path}")

        if versions_path.is_dir():
            version_files = [
                vp for vp in versions_path.iterdir() if vp.is_file()
            ]
        elif versions_path.is_file():
            version_files = [versions_path]
        else:
            raise IOError(f"{versions_path} should be either a directory or a file")
        for fname in version_files:
            versions = {}
            with open(fname, "r", encoding="utf-8") as fh:
                for line in fh.readlines():
                    name, version = line.strip().split(",")
                    versions[name] = version
        self.versions = versions
        return self.versions

    def to_json(self, filename):
        """Save class as JSON."""
        with open(filename, 'w') as f:
            json.dump(asdict(self), f, default=str, indent=2)
