import io
import matplotlib.pyplot as plt
import base64
import logomaker
import pandas as pd
import numpy as np
import warnings

plt.switch_backend('Agg')
def plot_motif(motif):
    motif_info = [
    motif.name,
    motif.matrix_id,
    motif.collection,
    motif.tf_class,
    motif.tf_family,
    motif.data_type,
    motif.medline
    ]
    df_counts = pd.DataFrame(motif.counts)
    # pseudocounts
    df_counts += 1
    # Calculate total counts per position
    total_counts = df_counts.sum(axis=1)
    # Calculate p_ai (normalized frequency per position)
    p_ai = df_counts.div(total_counts, axis=0)
    # Assume uniform background probabilities q_a
    q_a = 0.25

    log_ratios = np.log2(p_ai / q_a)
    E_Si = (p_ai * log_ratios).sum(axis=1)
    heights = p_ai.multiply(E_Si, axis=0)
    plt.figure(figsize=(10, 6))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*incompatible dtype.*")
        logo = logomaker.Logo(heights, shade_below=0.5, color_scheme="classic")
        
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    logo.ax.set_ylabel("Bits")
    logo.ax.set_ylim(0, 2)

    # Save the plot as a PNG image in memory
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png')
    img_buffer.seek(0)
    plt.close()

    img_data = base64.b64encode(img_buffer.read()).decode('utf-8')
    img_buffer.close()

    return motif_info,img_data