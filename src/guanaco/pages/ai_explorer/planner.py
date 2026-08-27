"""Question planning and strict validation for the AI explorer demo.

No expression or observation values are sent to the remote planner. Plotting is
performed locally after every returned field and plot type has been validated.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import requests
from pandas.api.types import is_bool_dtype, is_numeric_dtype


SUPPORTED_PLOTS = {"violin", "count_bar", "stacked_bar", "scatter", "heatmap"}
FIELD_KINDS = {"categorical", "numeric", "feature"}


@dataclass(frozen=True)
class Field:
    key: str
    label: str
    kind: str


@dataclass(frozen=True)
class DataProfile:
    fields: tuple[Field, ...]
    default_features: tuple[str, ...]
    n_observations: int
    n_features: int

    @property
    def by_key(self):
        return {field.key: field for field in self.fields}

    def for_prompt(self):
        return {
            "n_observations": self.n_observations,
            "n_features": self.n_features,
            "available_fields": [
                {"key": field.key, "label": field.label, "kind": field.kind}
                for field in self.fields
            ],
            "privacy_note": "No matrix or observation values are included.",
        }


def _categorical_obs(series: pd.Series) -> bool:
    if is_bool_dtype(series.dtype) or isinstance(series.dtype, pd.CategoricalDtype):
        return True
    if not is_numeric_dtype(series.dtype):
        return 1 < series.nunique(dropna=True) <= 40
    return False


def _question_feature_matches(adata, question: str) -> list[str]:
    text = question.casefold()
    return [
        str(name)
        for name in adata.var_names
        if _term_position(text, str(name)) is not None
    ]


def _term_position(text: str, term: str) -> int | None:
    """Find a complete field name, never a substring of another identifier."""
    match = re.search(
        rf"(?<![\w]){re.escape(term.casefold())}(?![\w])",
        text,
    )
    return match.start() if match else None


def build_data_profile(
    adata,
    question: str,
    default_features=None,
    *,
    feature_limit: int = 40,
) -> DataProfile:
    """Build a bounded schema profile without reading the expression matrix."""
    fields: list[Field] = []
    for name in adata.obs.columns:
        series = adata.obs[name]
        if _categorical_obs(series):
            fields.append(Field(f"obs:{name}", str(name), "categorical"))
        elif is_numeric_dtype(series.dtype):
            fields.append(Field(f"obs:{name}", str(name), "numeric"))

    requested = _question_feature_matches(adata, question)
    candidates = []
    valid_names = {str(name) for name in adata.var_names}
    for name in [*(default_features or ()), *requested]:
        name = str(name)
        if name in valid_names and name not in candidates:
            candidates.append(name)
    if not candidates:
        candidates.extend(str(name) for name in adata.var_names[: min(8, adata.n_vars)])
    candidates = candidates[:feature_limit]
    fields.extend(Field(f"feature:{name}", name, "feature") for name in candidates)
    return DataProfile(
        fields=tuple(fields),
        default_features=tuple(
            name for name in (default_features or ()) if str(name) in valid_names
        ),
        n_observations=int(adata.n_obs),
        n_features=int(adata.n_vars),
    )


PLAN_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["answerable", "answer_summary", "missing_data", "plots"],
    "properties": {
        "answerable": {"type": "boolean"},
        "answer_summary": {"type": "string"},
        "missing_data": {"type": "array", "items": {"type": "string"}},
        "plots": {
            "type": "array",
            "maxItems": 3,
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": [
                    "type",
                    "title",
                    "x",
                    "y",
                    "color",
                    "features",
                    "group_by",
                    "reason",
                ],
                "properties": {
                    "type": {"type": "string", "enum": sorted(SUPPORTED_PLOTS)},
                    "title": {"type": "string"},
                    "x": {"type": ["string", "null"]},
                    "y": {"type": ["string", "null"]},
                    "color": {"type": ["string", "null"]},
                    "features": {
                        "type": "array",
                        "items": {"type": "string"},
                        "maxItems": 10,
                    },
                    "group_by": {"type": ["string", "null"]},
                    "reason": {"type": "string"},
                },
            },
        },
    },
}


def _extract_response_text(response: dict[str, Any]) -> str:
    for item in response.get("output", []):
        for content in item.get("content", []):
            if content.get("type") == "output_text" and content.get("text"):
                return content["text"]
    raise ValueError("The LLM response did not contain structured output.")


def _remote_plan(question: str, profile: DataProfile) -> dict[str, Any]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not configured")
    instructions = (
        "You plan biological data visualizations. Use only fields listed in the "
        "provided data profile and copy their qualified keys exactly. First decide "
        "whether the data can answer the question. If a required measurement is "
        "absent, set answerable=false, name the missing data, and return no plots. "
        "Never infer that an unlisted field exists. Each plot must materially help "
        "answer the question. Features are numeric measurements. For violin use a "
        "categorical x, numeric/feature y, and optional categorical color; for "
        "stacked_bar use categorical x and color; for scatter use numeric/feature x "
        "and y with optional categorical color; for count_bar use one categorical x; "
        "for heatmap use feature keys and a categorical group_by. Do not claim "
        "causality or biological validation."
    )
    response = requests.post(
        "https://api.openai.com/v1/responses",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        json={
            "model": os.environ.get("GUANACO_LLM_MODEL", "gpt-5.4-mini"),
            "store": False,
            "input": [
                {"role": "developer", "content": instructions},
                {
                    "role": "user",
                    "content": json.dumps(
                        {"question": question, "data_profile": profile.for_prompt()}
                    ),
                },
            ],
            "text": {
                "format": {
                    "type": "json_schema",
                    "name": "guanaco_visualization_plan",
                    "strict": True,
                    "schema": PLAN_SCHEMA,
                }
            },
        },
        timeout=45,
    )
    response.raise_for_status()
    return json.loads(_extract_response_text(response.json()))


def _mentioned_fields(question: str, profile: DataProfile) -> list[Field]:
    text = question.casefold()
    matches = []
    for field in profile.fields:
        label = field.label.casefold()
        position = _term_position(text, label)
        if position is not None:
            matches.append((position, -len(label), field))
    # Respect question order, preferring the longest field at the same position
    # (e.g. IGHV_status before an IGHV feature contained in that name).
    matches.sort(key=lambda item: (item[0], item[1]))
    return [field for _position, _length, field in matches]


def local_demo_plan(question: str, profile: DataProfile) -> dict[str, Any]:
    """Conservative no-key fallback used only to make the demo immediately testable."""
    matched = _mentioned_fields(question, profile)
    categorical = [field for field in matched if field.kind == "categorical"]
    numeric = [field for field in matched if field.kind in {"numeric", "feature"}]
    lower = question.casefold()
    generic_feature_intent = any(
        term in lower
        for term in ("expression", "gene", "marker", "feature", "profile", "heatmap")
    )
    if (
        not numeric
        and categorical
        and profile.default_features
        and generic_feature_intent
    ):
        defaults = set(profile.default_features)
        numeric = [
            field
            for field in profile.fields
            if field.kind == "feature" and field.label in defaults
        ][:6]

    composition_intent = any(
        term in lower
        for term in ("composition", "proportion", "frequency", "count", "association", "relationship")
    )
    count_intent = any(
        term in lower
        for term in ("how many", "number of", "count", "frequency", "distribution")
    )
    plot = None
    if "heatmap" in lower and categorical and numeric:
        plot = {
            "type": "heatmap",
            "title": f"Group-level profile by {categorical[0].label}",
            "x": None,
            "y": None,
            "color": None,
            "features": [field.key for field in numeric[:8]],
            "group_by": categorical[0].key,
            "reason": "Compare the selected measurements across groups.",
        }
    elif numeric and categorical:
        plot = {
            "type": "violin",
            "title": f"{numeric[0].label} by {categorical[0].label}",
            "x": categorical[0].key,
            "y": numeric[0].key,
            "color": categorical[1].key if len(categorical) > 1 else None,
            "features": [],
            "group_by": None,
            "reason": "Show the distribution, group differences, and individual observations.",
        }
    elif categorical and count_intent:
        plot = {
            "type": "count_bar",
            "title": f"Number of observations by {categorical[0].label}",
            "x": categorical[0].key,
            "y": None,
            "color": None,
            "features": [],
            "group_by": None,
            "reason": "Count the observations in every category within this modality.",
        }
    elif len(categorical) >= 2 and composition_intent:
        plot = {
            "type": "stacked_bar",
            "title": f"{categorical[1].label} composition by {categorical[0].label}",
            "x": categorical[0].key,
            "y": None,
            "color": categorical[1].key,
            "features": [],
            "group_by": None,
            "reason": "Compare category composition between groups.",
        }
    elif len(numeric) >= 2:
        plot = {
            "type": "scatter",
            "title": f"{numeric[0].label} versus {numeric[1].label}",
            "x": numeric[0].key,
            "y": numeric[1].key,
            "color": None,
            "features": [],
            "group_by": None,
            "reason": "Inspect association and outlying observations.",
        }

    if plot is None:
        available = ", ".join(field.label for field in profile.fields[:8])
        return {
            "answerable": False,
            "answer_summary": (
                "I cannot map this question to the measurements in this data. "
                "Name at least one measured variable and, for a comparison, a group."
            ),
            "missing_data": [
                f"No requested measured variable was found. Examples: {available}."
            ],
            "plots": [],
        }
    return {
        "answerable": True,
        "answer_summary": "The required measurements are present in this modality.",
        "missing_data": [],
        "plots": [plot],
    }


def validate_plan(plan: dict[str, Any], profile: DataProfile) -> dict[str, Any]:
    """Reject hallucinated fields and semantically invalid plot specifications."""
    if not isinstance(plan, dict):
        raise ValueError("The planner did not return an object.")
    answerable = plan.get("answerable") is True
    plots = plan.get("plots")
    if not isinstance(plots, list):
        raise ValueError("The planner did not return a plot list.")
    if not answerable:
        return {
            "answerable": False,
            "answer_summary": str(plan.get("answer_summary") or "The data cannot answer this question."),
            "missing_data": [str(item) for item in plan.get("missing_data", [])],
            "plots": [],
        }
    if not plots:
        raise ValueError("The planner marked the question answerable but returned no plot.")

    fields = profile.by_key
    clean_plots = []
    for spec in plots[:3]:
        plot_type = spec.get("type")
        if plot_type not in SUPPORTED_PLOTS:
            raise ValueError(f"Unsupported plot type: {plot_type!r}")
        used = [spec.get("x"), spec.get("y"), spec.get("color"), spec.get("group_by")]
        used.extend(spec.get("features") or [])
        unknown = [key for key in used if key is not None and key not in fields]
        if unknown:
            raise ValueError(f"Planner requested unavailable fields: {', '.join(unknown)}")

        def kind(key):
            return fields[key].kind if key else None

        if plot_type == "violin" and not (
            kind(spec.get("x")) == "categorical"
            and kind(spec.get("y")) in {"numeric", "feature"}
            and kind(spec.get("color")) in {None, "categorical"}
        ):
            raise ValueError("A violin plot requires categorical x and numeric y.")
        if plot_type == "stacked_bar" and not (
            kind(spec.get("x")) == kind(spec.get("color")) == "categorical"
        ):
            raise ValueError("A stacked bar requires two categorical fields.")
        if plot_type == "count_bar" and not (
            kind(spec.get("x")) == "categorical"
            and spec.get("y") is None
            and spec.get("color") is None
        ):
            raise ValueError("A count bar requires one categorical field.")
        if plot_type == "scatter" and not (
            kind(spec.get("x")) in {"numeric", "feature"}
            and kind(spec.get("y")) in {"numeric", "feature"}
            and kind(spec.get("color")) in {None, "categorical"}
        ):
            raise ValueError("A scatter plot requires numeric x and y.")
        if plot_type == "heatmap" and not (
            kind(spec.get("group_by")) == "categorical"
            and spec.get("features")
            and all(kind(key) == "feature" for key in spec["features"])
        ):
            raise ValueError("A heatmap requires features and a categorical group.")
        clean_plots.append(dict(spec))
    return {
        "answerable": True,
        "answer_summary": str(plan.get("answer_summary") or "The data can address this question."),
        "missing_data": [],
        "plots": clean_plots,
    }


def plan_question(question: str, profile: DataProfile) -> tuple[dict[str, Any], str]:
    """Return a validated plan and the planner mode shown in the UI."""
    if os.environ.get("OPENAI_API_KEY"):
        try:
            return validate_plan(_remote_plan(question, profile), profile), "OpenAI planner"
        except Exception as exc:
            plan = validate_plan(local_demo_plan(question, profile), profile)
            return plan, f"Local demo planner (LLM unavailable: {type(exc).__name__})"
    return validate_plan(local_demo_plan(question, profile), profile), "Local demo planner"


def field_values(adata, field_key: str) -> pd.Series:
    """Resolve one already-validated field locally."""
    namespace, name = field_key.split(":", 1)
    if namespace == "obs":
        return pd.Series(adata.obs[name].to_numpy(), index=adata.obs_names, name=name)
    position = adata.var_names.get_loc(name)
    values = adata[:, position].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    return pd.Series(np.asarray(values).reshape(-1), index=adata.obs_names, name=name)
