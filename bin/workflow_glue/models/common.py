"""Common model classes used across all workflows."""
from dataclasses import asdict, dataclass, field
from decimal import Decimal
from enum import Enum
import json
from typing import Any, Dict, List


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
    barcode: str = field(
        metadata={
            "title": "Sample barcode",
            "description": "The physical barcode assigned to the sample"})
    sample_type: SampleType = field(
        metadata={
            "title": "Sample type",
            "description": "The type of the sample"})
    sample_pass: bool = field(
        metadata={
            "title": "Sample pass",
            "description": "If true the sample has passed workflow checks"})
    additional_identifiers: List[SampleIdentifier] = field(
        default_factory=list, metadata={
            "title": "Additional sample identifiers",
            "description": "Addition identifiers for the sample"})
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
        for indentifier in self.additional_identifiers:
            if indentifier.name == sample_identifier:
                return indentifier.value
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
            json.dump(asdict(self), f, default=str, indent=2, cls=DecimalEncoder)


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
class WorkflowResult():
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

    def to_json(self, filename):
        """Save class as JSON."""
        with open(filename, 'w') as f:
            json.dump(asdict(self), f, default=str, indent=2, cls=DecimalEncoder)


class DecimalEncoder(json.JSONEncoder):
    """This should probably be moved."""

    def default(self, obj):
        """Override the default method to handle Decimal objects."""
        if isinstance(obj, Decimal):
            return float(obj)
        return super().default(obj)
