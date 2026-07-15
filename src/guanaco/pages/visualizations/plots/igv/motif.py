"""Motif rendering helpers used by the IGV side panel."""

import base64
import io
import warnings

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend("Agg")

MOTIF_INFO_LABELS = [
    "TF Name",
    "Matrix ID",
    "Collection",
    "TF Class",
    "TF Family",
    "Data Type",
    "Medline",
]


def render_motif(motif):
    """Return ordered motif metadata and a base64-encoded sequence logo."""
    motif_info = [
        motif.name,
        motif.matrix_id,
        motif.collection,
        motif.tf_class,
        motif.tf_family,
        motif.data_type,
        motif.medline,
    ]
    counts = pd.DataFrame(motif.counts) + 1
    frequencies = counts.div(counts.sum(axis=1), axis=0)
    information = (frequencies * np.log2(frequencies / 0.25)).sum(axis=1)
    heights = frequencies.multiply(information, axis=0)

    plt.figure(figsize=(10, 6))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*incompatible dtype.*")
        logo = logomaker.Logo(
            heights,
            shade_below=0.5,
            color_scheme="classic",
        )
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.style_xticks(rotation=90, fmt="%d", anchor=0)
    logo.ax.set_ylabel("Bits")
    logo.ax.set_ylim(0, 2)

    image = io.BytesIO()
    plt.savefig(image, format="png")
    image.seek(0)
    plt.close()
    encoded = base64.b64encode(image.read()).decode("utf-8")
    image.close()
    return motif_info, encoded

