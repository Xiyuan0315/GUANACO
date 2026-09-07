"""Five reusable examples: one bundled public spatial dataset and four simulations."""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import anndata as ad
import guanaco as gc
import numpy as np
import pandas as pd
from scipy import sparse

from demo_data import (
    load_spatial_relationship_demo,
    make_mudata_demo,
    make_pathway_demo,
    make_single_cell,
)


def cell_summaries(adata=None):
    """Notebook demo 2: selected cells drive two expression summaries."""
    if adata is None:
        adata = make_single_cell(seed=202)
    markers = ["CD3D", "IL7R", "CCL5", "MS4A1", "NKG7", "LST1"]
    return gc.pl.linked_view(
        adata, prefix="gallery-cells", title="Selected cells → expression summaries",
        layout="grid",
        views=[
            gc.pl.umap(id="cells", color="cell_type", title="Select a cell population", height="400px"),
            gc.pl.dotplot(id="dotplot", var_names=markers, groupby="cell_type",
                          title="Expression fraction and mean", height="400px"),
            gc.pl.heatmap(id="heatmap", var_names=markers, groupby="cell_type",
                          title="Expression in selected cells", height="400px"),
        ],
        links=[gc.pl.link("cells", target, action="filter") for target in ("dotplot", "heatmap")],
    )


def rna_protein(mdata=None, rna_feature="CD4", protein_feature="CD4 protein"):
    """Notebook demo 4: shared cells connect different molecular measurements."""
    if mdata is None:
        mdata, rna_feature, protein_feature = make_mudata_demo()
    return gc.pl.linked_view(
        mdata, prefix="gallery-multiomics", title="RNA → protein measurements",
        views=[
            gc.pl.umap(id="rna", data="rna", color=rna_feature,
                       title=f"RNA · {rna_feature}", height="440px"),
            gc.pl.pca(id="protein", data="protein", color=protein_feature,
                      title=f"Protein PCA · {protein_feature}", height="440px"),
        ],
        links=[gc.pl.link("rna", "protein")],
    )


def spatial_relationships(data=None):
    """Notebook demo 8: one relationship maps to both cells and distance records."""
    if data is None:
        path = Path(__file__).with_name("data") / "visium_hne_spatial.h5ad"
        if not path.is_file():
            raise FileNotFoundError(
                f"Bundled spatial data missing: {path}. Run prepare_spatial.py with "
                "the original notebook dataset; synthetic fallback is not used here."
            )
        pairs, spatial, curves, _ = load_spatial_relationship_demo(path)
        data = {"pairs": pairs, "spatial": spatial, "cooccurrence": curves}
    groups = list(data["spatial"].obs["cluster"].cat.categories)
    return gc.pl.linked_view(
        data, prefix="gallery-spatial", title="Spatial relationships → locations and distance",
        layout="grid",
        views=[
            gc.pl.view("plotly.heatmap", id="neighborhoods", data="pairs",
                       x="target_group", y="source_group", value="enrichment",
                       color_map="RdBu_r", colorbar_title="Z-score",
                       x_order=groups, y_order=groups,
                       x_title="Neighbor group", y_title="Conditional group",
                       title="Original neighborhood enrichment · click a pair", height="560px"),
            gc.pl.spatial(id="locations", data="spatial", color="cluster", size=3,
                          title="Original H&E image and spatial spots", height="560px"),
            gc.pl.view("plotly.line", id="cooccurrence", data="cooccurrence",
                       x="distance", y="co_occurrence", group="pair_label",
                       color="pair_label", color_mode="categorical", fixed_axes=True,
                       title="Co-occurrence across distance", height="430px"),
        ],
        links=[
            gc.pl.link("neighborhoods", "locations", by="cell", key="cell_id", action="filter"),
            gc.pl.link("neighborhoods", "cooccurrence", key="pair_id"),
        ],
    )


def pathway_genes(data=None):
    """Notebook demo 10: a pathway selection expands to its constituent genes."""
    if data is None:
        pathways, genes = make_pathway_demo()
        data = {"pathways": pathways, "genes": genes}
    return gc.pl.linked_view(
        data, prefix="gallery-pathways", title="Pathways → gene expression",
        views=[
            gc.pl.view("plotly.bar", id="pathways", data="pathways",
                       x="pathway_score", y="pathway", orientation="h", color="pathway",
                       title="Click a pathway", height="420px"),
            gc.pl.view("plotly.heatmap", id="genes", data="genes",
                       x="gene", y="cell_type", value="gene_mean", color_map="Blues",
                       colorbar_title="Mean expression", title="Genes in the selected pathway",
                       height="460px"),
        ],
        links=[gc.pl.link("pathways", "genes", key="pathway")],
    )


def make_peak_data():
    """Invented coordinates and counts, NOT a human genome reference or experiment."""
    rna = make_single_cell(360, seed=111)
    genes = ["CD4", "IL7R", "MS4A1", "NKG7", "LST1"]
    rng = np.random.default_rng(112)
    starts = np.concatenate([100_000 + i * 200_000 + np.arange(18) * 1_000 for i in range(5)])
    activity = np.asarray(rna[:, genes].X)
    probabilities = np.repeat(0.04 + 0.7 * activity / (activity.max(axis=0) + 0.1), 18, axis=1)
    counts = (rng.random(probabilities.shape) < probabilities).astype(np.float32)
    atac = ad.AnnData(
        X=sparse.csr_matrix(counts), obs=rna.obs.copy(),
        var=pd.DataFrame(index=[f"chrDemo:{start}-{start + 600}" for start in starts]),
    )
    annotation = Path(__file__).with_name("synthetic_genes.gtf")
    return {"rna": rna, "atac": atac}, genes, annotation


def gene_accessibility(data=None, genes=None, annotation=None):
    """Notebook demo 11: annotation maps a gene name to an ATAC region."""
    if data is None:
        data, genes, annotation = make_peak_data()
    return gc.pl.linked_view(
        data, prefix="gallery-atac", title="Gene expression → chromatin accessibility",
        views=[
            gc.pl.matrixplot(id="gene_overview", data="rna", var_names=genes,
                             groupby="cell_type", standardization="var",
                             title="Grouped expression · click a gene", height="560px"),
            gc.pl.peak_browser(id="peaks", data="atac", region=genes[0],
                               gene_annotation=str(annotation), groupby="cell_type",
                               metric="detection", y_mode="shared",
                               title="Accessibility around the selected gene", height="620px"),
        ],
        links=[gc.pl.link("gene_overview", "peaks", by="feature")],
    )


@dataclass(frozen=True)
class Case:
    slug: str
    title: str
    question: str
    instruction: str
    explanation: str
    data_note: str
    build: Callable


CASES = (
    Case("cells", "Cell populations",
         "How does a selected population differ in its molecular profile?",
         "Draw a lasso around cells in the embedding. Both expression views update to your selection.",
         "One cell selection filters two complementary summaries: expression prevalence and individual-cell expression.",
         "540 simulated cells, six cell types and 16 genes. The embedding is simulated, not a computed UMAP.", cell_summaries),
    Case("rna-protein", "RNA & protein",
         "Where do RNA-selected cells appear in protein space?",
         "Draw a lasso around cells in the RNA plot. Matching cells are highlighted in the protein PCA.",
         "Shared cell identifiers connect modalities even when their features and coordinates differ.",
         "420 simulated cells with paired RNA and protein measurements; simulated RNA embedding and computed protein PCA.", rna_protein),
    Case("spatial", "Spatial relationships",
         "Where are two cell groups located, and how does their association vary with distance?",
         "Click an off-diagonal tile, for example Hippocampus × Pyramidal_layer. The tissue overlay filters to those groups and the distance curve updates.",
         "The selected pair maps to cell identifiers for the spatial view and a pair identifier for the distance curve.",
         "Original Visium mouse-brain dataset from the Squidpy H&E tutorial (10x Genomics). All 2,688 spots, 15 annotated groups, original H&E images and scale factors, the complete 15 × 15 enrichment matrix, and precomputed co-occurrence curves are retained. Only unused expression data and metadata are omitted. Spatial spots can contain multiple cells.", spatial_relationships),
    Case("pathways", "Pathways & genes",
         "Which genes contribute to the expression pattern of a pathway?",
         "Click a pathway bar to show its genes across cell types.",
         "A shared pathway key connects one overview mark to multiple gene records.",
         "Illustrative pathway memberships and mean-expression scores from simulated cells; not enrichment-test results.", pathway_genes),
    Case("atac", "Genes & accessibility",
         "What does chromatin accessibility look like around a selected gene?",
         "Click a gene tile. The accessibility browser moves to that gene’s region.",
         "Gene annotation translates the selected gene into genomic coordinates; RNA and ATAC do not need identical feature names.",
         "360 simulated cells and 90 peaks on an invented chromosome (chrDemo). Gene coordinates and counts are synthetic.", gene_accessibility),
)
