import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")                        
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO
from hotspot_finder import HotspotRegion
COLOURS = {
    "POTENTIAL PCR ARTEFACT":           "#E85D24",   
    "HIGHLY POTENTIAL NATURAL MUTATION": "#1D9E75",  
    "LOW POTENTIAL NATURAL MUTATION":    "#F2A623",  
    "HIGH":                             "#A32D2D",   
    "MODERATE":                         "#F2A623",   
    "homopolymer":                      "#3B8BD4",   
    "STR":                              "#7F77DD",   
}
def generate_all_reports(
    classified_df: pd.DataFrame,
    hotspot_map:   Dict[str, List[HotspotRegion]],
    reference:     str,
    output_dir:    str,
    logger         = None
):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    if logger:
        logger.info(f"  Output directory: {out}")
    genome_sizes = {rec.id: len(rec.seq)
                    for rec in SeqIO.parse(reference, "fasta")}
    full_tsv = out / "variants_classified.tsv"
    if not classified_df.empty:
        classified_df.to_csv(full_tsv, sep="\t", index=False)
        if logger: logger.info(f"  Written: {full_tsv.name}")
    _write_per_site_counts(classified_df, out / "per_site_counts.tsv", logger)
    artefacts = classified_df[
        classified_df["classification"] == "POTENTIAL PCR ARTEFACT"
    ] if not classified_df.empty else pd.DataFrame()
    artefacts.to_csv(out / "artefacts_flagged.tsv", sep="\t", index=False)
    if logger: logger.info(f"  Written: artefacts_flagged.tsv  ({len(artefacts)} sites)")
    naturals = classified_df[
        classified_df["classification"].str.contains("NATURAL")
    ] if not classified_df.empty else pd.DataFrame()
    naturals.to_csv(out / "natural_mutations.tsv", sep="\t", index=False)
    if logger: logger.info(f"  Written: natural_mutations.tsv  ({len(naturals)} sites)")
    _write_hotspot_summary(hotspot_map, out / "hotspot_map_summary.tsv", logger)
    _write_text_summary(classified_df, hotspot_map, out / "summary_report.txt", logger)
    _plot_hotspot_positions(hotspot_map, genome_sizes,
                            out / "figures" / "hotspot_positions.png")
    if not classified_df.empty:
        _plot_classification_overview(classified_df, out / "figures" / "classification_overview.png")
        _plot_quality_distribution(classified_df, out / "figures" / "quality_distribution.png")
        _plot_genome_map(classified_df, hotspot_map, genome_sizes,
                         out / "figures" / "genome_map.png")
        if "AG_Severity_Score" in classified_df.columns and not classified_df["AG_Severity_Score"].isna().all():
            _plot_alphagenome_impact(classified_df, out / "figures" / "alphagenome_impact.png")
        if logger: logger.info("  Figures written to figures/")
    if logger:
        logger.info(f"  All reports complete.")
def _write_per_site_counts(df: pd.DataFrame, path: Path, logger):
    if df.empty:
        pd.DataFrame(columns=[
            "chrom","pos_1","n_insertions","n_deletions","total_indels",
            "classification","mean_vaf","mean_quality"
        ]).to_csv(path, sep="\t", index=False)
        return
    grouped = (
        df.groupby(["chrom", "pos_1", "classification"])
        .agg(
            n_insertions = ("mut_type", lambda x: (x == "insertion").sum()),
            n_deletions  = ("mut_type", lambda x: (x == "deletion").sum()),
            total_indels = ("mut_type", "count"),
            mean_vaf     = ("vaf_percent", "mean"),
            mean_quality = ("mean_alt_supporting_phred", "mean"),
        )
        .reset_index()
        .sort_values(["chrom", "pos_1"])
    )
    grouped["mean_vaf"]     = grouped["mean_vaf"].round(2)
    grouped["mean_quality"] = grouped["mean_quality"].round(1)
    grouped.to_csv(path, sep="\t", index=False)
    if logger: logger.info(f"  Written: {path.name}  ({len(grouped)} positions)")
def _write_hotspot_summary(hotspot_map, path, logger):
    rows = []
    for contig, regions in hotspot_map.items():
        for r in regions:
            rows.append({
                "contig":            r.contig,
                "start_1based":      r.start + 1,
                "end_1based":        r.end,
                "total_length_bp":   r.total_length,
                "region_type":       r.region_type,
                "repeat_unit":       r.unit,
                "unit_length_bp":    r.unit_length,
                "copies":            round(r.copies, 1),
                "risk_level":        r.risk_level,
                "description":       r.description,
            })
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    if logger: logger.info(f"  Written: {path.name}  ({len(df)} hotspot regions)")
def _write_text_summary(df, hotspot_map, path, logger):
    n_total      = len(df)
    n_artefact   = (df["classification"] == "POTENTIAL PCR ARTEFACT").sum() if n_total > 0 else 0
    n_high_nat   = (df["classification"] == "HIGHLY POTENTIAL NATURAL MUTATION").sum() if n_total > 0 else 0
    n_low_nat    = (df["classification"] == "LOW POTENTIAL NATURAL MUTATION").sum()    if n_total > 0 else 0
    n_insertions = (df["mut_type"] == "insertion").sum() if n_total > 0 else 0
    n_deletions  = (df["mut_type"] == "deletion").sum()  if n_total > 0 else 0
    total_hotspots = sum(len(v) for v in hotspot_map.values())
    high_hotspots  = sum(1 for v in hotspot_map.values() for r in v if r.risk_level == "HIGH")
    with open(path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("  PCR SLIPPAGE ARTEFACT DETECTION — SUMMARY REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write("REFERENCE SLIPPAGE HOTSPOT MAP\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Total hotspot regions   : {total_hotspots}\n")
        f.write(f"  HIGH-risk regions       : {high_hotspots}\n")
        f.write(f"  MODERATE-risk regions   : {total_hotspots - high_hotspots}\n\n")
        f.write("INDEL VARIANT COUNTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Total indels called     : {n_total}\n")
        f.write(f"  Insertions              : {n_insertions}\n")
        f.write(f"  Deletions               : {n_deletions}\n\n")
        f.write("CLASSIFICATION RESULTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"  POTENTIAL PCR ARTEFACT             : {n_artefact}\n")
        f.write(f"  HIGHLY POTENTIAL NATURAL MUTATION  : {n_high_nat}\n")
        f.write(f"  LOW POTENTIAL NATURAL MUTATION     : {n_low_nat}\n\n")
        f.write("CLASSIFICATION DEFINITIONS\n")
        f.write("-" * 40 + "\n")
        f.write(
            "  POTENTIAL PCR ARTEFACT\n"
            "    The indel falls within a homopolymer or STR region of the reference.\n"
            "    Polymerase slippage during PCR amplification is a known cause of\n"
            "    false indels at such positions. These should NOT be treated as\n"
            "    real mutations without additional experimental validation.\n\n"
        )
        f.write(
            "  HIGHLY POTENTIAL NATURAL MUTATION\n"
            "    The indel is outside any slippage-prone region AND is supported\n"
            "    by reads with high mean Phred quality scores. This combination\n"
            "    strongly suggests the indel is a genuine biological event.\n\n"
        )
        f.write(
            "  LOW POTENTIAL NATURAL MUTATION\n"
            "    The indel is outside any slippage-prone region BUT is supported\n"
            "    by reads with low Phred quality scores. The indel may be real\n"
            "    but should be treated with caution — could be a sequencing error.\n\n"
        )
        if n_total > 0:
            f.write("TOP FLAGGED SITES (potential PCR artefacts at highest-risk hotspots)\n")
            f.write("-" * 40 + "\n")
            artefact_df = df[df["classification"] == "POTENTIAL PCR ARTEFACT"]
            if not artefact_df.empty:
                top = artefact_df.nlargest(10, "total_depth")
                for _, row in top.iterrows():
                    f.write(
                        f"  {row['chrom']}:{row['pos_1']:<8} "
                        f"{row['mut_type'].upper()} {row['mut_size_bp']}bp  "
                        f"VAF={row['vaf_percent']:.1f}%  "
                        f"Depth={row['total_depth']}  "
                        f"Hotspot=({row['hotspot_unit']})x{row['hotspot_copies']:.0f} "
                        f"[{row['hotspot_risk_level']}]\n"
                    )
            f.write("\n")
            f.write("TOP NATURAL MUTATION CANDIDATES\n")
            f.write("-" * 40 + "\n")
            nat_df = df[df["classification"] == "HIGHLY POTENTIAL NATURAL MUTATION"]
            if not nat_df.empty:
                top_nat = nat_df.nlargest(10, "mean_alt_supporting_phred")
                for _, row in top_nat.iterrows():
                    f.write(
                        f"  {row['chrom']}:{row['pos_1']:<8} "
                        f"{row['mut_type'].upper()} {row['mut_size_bp']}bp  "
                        f"VAF={row['vaf_percent']:.1f}%  "
                        f"Depth={row['total_depth']}  "
                        f"MeanPhred={row['mean_alt_supporting_phred']:.1f}\n"
                    )
    if logger: logger.info(f"  Written: {path.name}")
def _plot_classification_overview(df, path):
    fig, (ax_pie, ax_bar) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Indel Classification Overview", fontsize=13, fontweight="bold")
    counts = df["classification"].value_counts()
    labels = list(counts.index)
    sizes  = list(counts.values)
    colors = [COLOURS.get(l, "#888") for l in labels]
    wedge_props = dict(width=0.5, edgecolor="white", linewidth=1.5)
    ax_pie.pie(sizes, labels=None, colors=colors, autopct="%1.1f%%",
               wedgeprops=wedge_props, startangle=140, pctdistance=0.75)
    ax_pie.set_title("Proportion", fontsize=11)
    patches = [mpatches.Patch(color=c, label=l) for l, c in zip(labels, colors)]
    ax_pie.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5, -0.2),
                  fontsize=8, frameon=False)
    ins_counts = df[df["mut_type"] == "insertion"]["classification"].value_counts()
    del_counts = df[df["mut_type"] == "deletion"]["classification"].value_counts()
    all_labels = sorted(set(list(ins_counts.index) + list(del_counts.index)))
    x = np.arange(len(all_labels))
    w = 0.35
    ins_vals = [ins_counts.get(l, 0) for l in all_labels]
    del_vals = [del_counts.get(l, 0) for l in all_labels]
    ax_bar.bar(x - w/2, ins_vals, w, label="Insertions", color="#3B8BD4", alpha=0.85)
    ax_bar.bar(x + w/2, del_vals, w, label="Deletions",  color="#D85A30", alpha=0.85)
    ax_bar.set_xticks(x)
    short_labels = [l.replace("POTENTIAL ", "").replace(" MUTATION", "")
                    .replace("PCR ARTEFACT", "PCR\nARTEFACT") for l in all_labels]
    ax_bar.set_xticklabels(short_labels, fontsize=8)
    ax_bar.set_ylabel("Count")
    ax_bar.set_title("Insertions vs Deletions by Class", fontsize=11)
    ax_bar.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
def _plot_quality_distribution(df, path):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Supporting-Read Quality by Classification", fontsize=12, fontweight="bold")
    classes   = [
        "POTENTIAL PCR ARTEFACT",
        "HIGHLY POTENTIAL NATURAL MUTATION",
        "LOW POTENTIAL NATURAL MUTATION"
    ]
    short_cls = ["PCR Artefact", "High Natural", "Low Natural"]
    data      = [df[df["classification"] == c]["mean_alt_supporting_phred"].dropna().values
                 for c in classes]
    colors    = [COLOURS[c] for c in classes]
    parts = ax.violinplot(
        [d if len(d) > 0 else [0] for d in data],
        positions=range(len(classes)),
        showmedians=True,
        showextrema=True
    )
    for pc, col in zip(parts["bodies"], colors):
        pc.set_facecolor(col)
        pc.set_alpha(0.6)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(2)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(short_cls, fontsize=10)
    ax.set_ylabel("Mean Phred Quality Score")
    ax.set_xlabel("Classification")
    ax.axhline(y=30, color="gray", linestyle="--", linewidth=0.8, label="Q30 threshold")
    ax.axhline(y=20, color="lightgray", linestyle="--", linewidth=0.8, label="Q20 threshold")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
def _plot_genome_map(df, hotspot_map, genome_sizes, path):
    contigs = list(genome_sizes.keys())
    n_rows  = len(contigs)
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 2.5 * n_rows + 1), squeeze=False)
    fig.suptitle("Genome-Wide Indel Map\n(background = slippage hotspots)",
                 fontsize=12, fontweight="bold")
    for row_idx, contig in enumerate(contigs):
        ax       = axes[row_idx][0]
        gen_len  = genome_sizes[contig]
        regions  = hotspot_map.get(contig, [])
        ctg_df   = df[df["chrom"] == contig] if not df.empty else pd.DataFrame()
        for reg in regions:
            band_col = "#FAEEDA" if reg.risk_level == "MODERATE" else "#FCEBEB"
            ax.axvspan(reg.start, reg.end, alpha=0.6, color=band_col, zorder=1)
        if not ctg_df.empty:
            for _, row in ctg_df.iterrows():
                col  = COLOURS.get(row["classification"], "#888")
                ax.vlines(row["pos_1"], 0, 1, colors=col, alpha=0.7,
                          linewidth=1.2, zorder=3)
        ax.set_xlim(0, gen_len)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Genomic Position (bp)", fontsize=8)
        ax.set_title(f"Contig: {contig}  ({gen_len:,} bp)", fontsize=9, loc="left")
    legend_elements = [
        mpatches.Patch(color=COLOURS["POTENTIAL PCR ARTEFACT"],          label="PCR Artefact"),
        mpatches.Patch(color=COLOURS["HIGHLY POTENTIAL NATURAL MUTATION"],label="High Natural"),
        mpatches.Patch(color=COLOURS["LOW POTENTIAL NATURAL MUTATION"],   label="Low Natural"),
        mpatches.Patch(color="#FAEEDA", label="MODERATE hotspot"),
        mpatches.Patch(color="#FCEBEB", label="HIGH hotspot"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=5,
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
def _plot_hotspot_positions(hotspot_map, genome_sizes, path):
    contigs = list(genome_sizes.keys())
    n_rows  = len(contigs)
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 2.8 * n_rows + 1), squeeze=False)
    fig.suptitle("Slippage Hotspot Density Across Reference Genome",
                 fontsize=12, fontweight="bold")
    window = 200   
    for row_idx, contig in enumerate(contigs):
        ax      = axes[row_idx][0]
        gen_len = genome_sizes[contig]
        regions = hotspot_map.get(contig, [])
        if not regions:
            ax.set_title(f"{contig} — no hotspots", fontsize=9)
            continue
        bins      = np.arange(0, gen_len + window, window)
        hp_counts = np.zeros(len(bins) - 1)
        st_counts = np.zeros(len(bins) - 1)
        for r in regions:
            bin_idx = min(r.start // window, len(bins) - 2)
            if r.region_type == "homopolymer":
                hp_counts[bin_idx] += 1
            else:
                st_counts[bin_idx] += 1
        bin_centres = (bins[:-1] + bins[1:]) / 2
        ax.bar(bin_centres, hp_counts, width=window * 0.85,
               color=COLOURS["homopolymer"], alpha=0.75, label="Homopolymer")
        ax.bar(bin_centres, st_counts, width=window * 0.85,
               bottom=hp_counts, color=COLOURS["STR"], alpha=0.75, label="STR")
        ax.set_xlim(0, gen_len)
        ax.set_ylabel("Count")
        ax.set_xlabel("Position (bp)", fontsize=8)
        ax.set_title(f"Contig: {contig}  (window = {window} bp)", fontsize=9, loc="left")
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
def _plot_alphagenome_impact(df, path):
    ag_df = df[df["AG_Severity_Score"].notna()]
    if ag_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle("AlphaGenome Functional Impact (Genuine Mutations)", fontsize=12, fontweight="bold")
    scatter = ax.scatter(
        ag_df["vaf_percent"], 
        ag_df["AG_Severity_Score"],
        c=ag_df["AG_Severity_Score"].astype(float),
        cmap="viridis",
        s=100,
        alpha=0.8,
        edgecolors="w",
        linewidth=0.5
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Predicted Severity Score")
    ax.set_xlabel("Variant Allele Frequency (%)")
    ax.set_ylabel("AlphaGenome Severity Score")
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=0.8, color="red", linestyle="--", linewidth=1, alpha=0.5, label="High Severity Threshold")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
