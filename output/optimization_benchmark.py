"""Compare current code with 447813f; synthetic data, cold caches, no network.

Run: .pixi/envs/default/bin/python output/optimization_benchmark.py
Timing is the median of five runs. Tracemalloc peak is measured separately and
excludes the already loaded input matrix, imports and index construction.
"""
import gc
import json
import subprocess
import sys
import time
import tracemalloc
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from plotly.io import to_json

from guanaco.pages.matrix.plots import atac_browser as atac, dotmatrix as dot
from guanaco.pages.matrix.plots.violin1 import plot_violin1
from guanaco.utils import gene_extraction_utils as genes


def baseline(name, path):
    module = types.ModuleType(name)
    sys.modules[name] = module
    code = subprocess.check_output(["git", "show", f"447813f:{path}"], text=True)
    exec(compile(code, f"447813f/{path}", "exec"), module.__dict__)
    return module


old_genes = baseline("baseline_genes", "src/guanaco/utils/gene_extraction_utils.py")
old_dot = baseline("baseline_dot", "src/guanaco/pages/matrix/plots/dotmatrix.py")
old_atac = baseline("baseline_atac", "src/guanaco/pages/matrix/plots/atac_browser.py")
old_dot.extract_gene_expression = old_genes.extract_gene_expression
old_dot.prewarm_gene_cache = old_genes.prewarm_gene_cache


def reset():
    genes.clear_gene_cache()
    old_genes.clear_gene_cache()
    for module in (atac, old_atac):
        for name in ("_signal_cache", "_col_cache", "_dense_cache"):
            cache = getattr(module, name, None)
            if cache is not None:
                cache._store.clear()
                if hasattr(cache, "bytes"):
                    cache.bytes = 0
    gc.collect()


def measure(label, fn):
    elapsed = []
    for _ in range(5):
        reset()
        start = time.perf_counter()
        result = fn()
        elapsed.append(time.perf_counter() - start)
        del result
    reset()
    tracemalloc.start()
    result = fn()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(json.dumps({"case": label, "median_ms": round(np.median(elapsed) * 1000, 2),
                      "peak_MiB": round(peak / 1024**2, 2)}), flush=True)
    return result


rng = np.random.default_rng(42)
n, p = 200_000, 80
x = sparse.random(n, p, density=0.01, format="csc", dtype=np.float32,
                  random_state=rng, data_rvs=lambda size: rng.integers(1, 10, size=size))
data = ad.AnnData(x, obs=pd.DataFrame({"group": pd.Categorical([f"G{i % 20}" for i in range(n)])},
                                    index=[f"c{i}" for i in range(n)]),
                  var=pd.DataFrame(index=[f"chr1:{i*100+1}-{i*100+80}" for i in range(p)]))
region = {"chrom": "chr1", "start": 0, "end": 8000}
for module in (atac, old_atac):
    module.build_peak_index(data)

# Warm plotting machinery, not the benchmark caches.
dot.plot_dot_matrix(data[:100].copy(), data.var_names[:2], "group", None)
for fmt in ("csc", "csr", "dense"):
    data.X = data.X.toarray() if fmt == "dense" else data.X.asformat(fmt)
    for name, before, after in (("dot", old_dot.plot_dot_matrix, dot.plot_dot_matrix),
                                 ("ATAC", old_atac.compute_atac_signal, atac.compute_atac_signal)):
        args = (data, data.var_names, "group", None) if name == "dot" else (data, region)
        kwargs = {} if name == "dot" else {"groupby": "group"}
        a = measure(f"{fmt}/{name}/before", lambda: before(*args, **kwargs))
        b = measure(f"{fmt}/{name}/after", lambda: after(*args, **kwargs))
        if name == "dot":
            np.testing.assert_allclose(a.data[0].marker.color, b.data[0].marker.color, rtol=2e-6)
            np.testing.assert_allclose(a.data[0].marker.size, b.data[0].marker.size, rtol=2e-6)
        else:
            assert [s["name"] for s in a["signals"]] == [s["name"] for s in b["signals"]]
            np.testing.assert_allclose([s["values"] for s in a["signals"]],
                                       [s["values"] for s in b["signals"]], rtol=2e-6)

# Small dense/sparse cases verify every existing transformation/standardization.
for fmt in ("dense", "csr"):
    small = data[:1000, :5].copy()
    if fmt == "dense":
        small.X = np.asarray(small.X)
    else:
        small.X = sparse.csr_matrix(small.X)
    for transform in (None, "log1p", "zscore"):
        for scale in (None, "zscore", "minmax", "group"):
            reset()
            kwargs = dict(transformation=transform, standardization=scale,
                          plot_type="matrixplot", selected_cells=small.obs_names[:753].tolist())
            a = old_dot.plot_dot_matrix(small, small.var_names, "group", ["G2", "G1"], **kwargs)
            b = dot.plot_dot_matrix(small, small.var_names, "group", ["G2", "G1"], **kwargs)
            np.testing.assert_allclose(a.data[0].z, b.data[0].z, rtol=2e-5, atol=1e-6)
print("Dense/sparse dot transformation and standardization equivalence: passed", flush=True)

# Consume >64MiB of requested columns without retaining the output, as heatmap does.
large = ad.AnnData(sparse.random(1_000_000, 24, density=0.01, format="csc", dtype=np.float32,
                               random_state=rng, data_rvs=lambda size: rng.integers(1, 10, size=size)))
for module, name in ((old_genes, "before"), (genes, "after")):
    def consume(module=module):
        if module is old_genes:
            module.prewarm_gene_cache(large, large.var_names.tolist())
            return sum(float(module.extract_gene_expression(large, g).sum()) for g in large.var_names)
        return sum(float(v.sum()) for _, v in module.iter_gene_expression(large, large.var_names))
    measure(f"1M/24genes/{name}", consume)
    reads = []
    vector, block = module._compute_gene_vector, module._compute_gene_block
    def counted_vector(adata, gene, *args):
        reads.append(1)
        return vector(adata, gene, *args)
    def counted_block(adata, names, *args, **kwargs):
        reads.append(len(names))
        return block(adata, names, *args, **kwargs)
    module._compute_gene_vector, module._compute_gene_block = counted_vector, counted_block
    reset()
    consume()
    print(json.dumps({"case": name, "read_calls": len(reads), "columns_read": sum(reads)}), flush=True)
    module._compute_gene_vector, module._compute_gene_block = vector, block

# Figure-bearing JSON only, before compression: prior cache -> returned cache ->
# display callback cache input -> figure output. Exclude old full-figure State,
# since that was already removed in the preceding optimization round.
violin_data = data[:10_000, :8].copy()
figure = json.loads(to_json(plot_violin1(violin_data, violin_data.var_names.tolist(), "group",
                                       labels=violin_data.obs["group"].cat.categories.tolist())))
encode = lambda value: len(json.dumps(value, separators=(",", ":")).encode())
previous = {str(i): figure for i in range(8)}
updated = {str(i): figure for i in range(9)}
before = encode(previous) + 2 * encode(updated) + encode(figure)
after = encode(figure)
print(json.dumps({"case": "violin/8to9history", "before_MiB": round(before / 1024**2, 2),
                  "after_MiB": round(after / 1024**2, 2), "reduction_pct": round(100*(1-after/before), 2)}))
