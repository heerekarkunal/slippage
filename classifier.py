import pysam
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from hotspot_finder import HotspotRegion, position_in_hotspot
LABEL_PCR_ARTEFACT  = "POTENTIAL PCR ARTEFACT"
LABEL_HIGH_NATURAL  = "HIGHLY POTENTIAL NATURAL MUTATION"
LABEL_LOW_NATURAL   = "LOW POTENTIAL NATURAL MUTATION"
def _parse_vcf(vcf_path: str, min_vaf: float = 0.01) -> List[dict]:
    import subprocess
    query_fmt = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t[%AD]\n"
    result = subprocess.run(
        ["bcftools", "query", "-f", query_fmt, vcf_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"bcftools query failed: {result.stderr}")
    variants = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, pos_str, ref, alt, dp_str = parts[0], parts[1], parts[2], parts[3], parts[4]
        ad_str = parts[5] if len(parts) > 5 else "."
        try:
            pos = int(pos_str) - 1  
            dp  = int(dp_str)  if dp_str not in (".", "") else 0
        except ValueError:
            continue
        ref_depth, alt_depth = 0, 0
        if ad_str not in (".", ""):
            ad_parts = ad_str.split(",")
            try:
                ref_depth = int(ad_parts[0])
                alt_depth = int(ad_parts[1]) if len(ad_parts) > 1 else 0
            except (ValueError, IndexError):
                pass
        total = ref_depth + alt_depth
        vaf   = alt_depth / total if total > 0 else 0.0
        if vaf < min_vaf:
            continue
        ref_len = len(ref)
        alt_len = len(alt)
        if alt_len > ref_len:
            mut_type = "insertion"
            mut_size = alt_len - ref_len
        elif alt_len < ref_len:
            mut_type = "deletion"
            mut_size = ref_len - alt_len
        else:
            continue  
        variants.append({
            "chrom":     chrom,
            "pos_0":     pos,          
            "pos_1":     pos + 1,      
            "ref":       ref,
            "alt":       alt,
            "depth":     dp,
            "ref_depth": ref_depth,
            "alt_depth": alt_depth,
            "vaf":       round(vaf, 4),
            "mut_type":  mut_type,
            "mut_size":  mut_size,
        })
    return variants
def _get_alt_supporting_quality(
    bam_path: str,
    chrom:    str,
    pos_0:    int,
    alt:      str,
    ref:      str,
    window:   int = 2
) -> Tuple[float, int]:
    qualities = []
    n_alt     = 0
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for pileup_col in bam.pileup(
            contig    = chrom,
            start     = max(0, pos_0 - window),
            stop      = pos_0 + len(ref) + window,
            truncate  = True,
            min_base_quality = 0,   
            stepper   = "all"
        ):
            if pileup_col.reference_pos != pos_0:
                continue
            for pileup_read in pileup_col.pileups:
                aln = pileup_read.alignment
                if aln.is_secondary or aln.is_supplementary:
                    continue
                qpos = pileup_read.query_position
                if qpos is None:
                    continue
                q_start = max(0, qpos - window)
                q_end   = min(len(aln.query_qualities or []), qpos + window + 1)
                if aln.query_qualities is not None:
                    window_quals = list(aln.query_qualities[q_start:q_end])
                    if window_quals:
                        qualities.extend(window_quals)
                        n_alt += 1
        bam.close()
    except Exception:
        pass  
    if qualities:
        return round(float(np.mean(qualities)), 1), n_alt
    return 0.0, n_alt
def classify_variants(
    vcf_path:              str,
    bam_path:              str,
    hotspot_map:           Dict[str, List[HotspotRegion]],
    high_quality_threshold: int = 30,
    min_quality:           int  = 20,
    min_vaf:               float = 0.01,
    logger                 = None
) -> pd.DataFrame:
    if logger:
        logger.info(f"  Reading variants from {Path(vcf_path).name}...")
    variants = _parse_vcf(vcf_path, min_vaf=min_vaf)
    if logger:
        logger.info(f"  Classifying {len(variants):,} indels...")
    rows = []
    for i, v in enumerate(variants):
        chrom = v["chrom"]
        pos_0 = v["pos_0"]
        hotspot = position_in_hotspot(chrom, pos_0, hotspot_map)
        if hotspot is None and v["mut_type"] == "deletion":
            for offset in range(1, v["mut_size"] + 1):
                hotspot = position_in_hotspot(chrom, pos_0 + offset, hotspot_map)
                if hotspot:
                    break
        mean_qual, n_alt_reads = _get_alt_supporting_quality(
            bam_path = bam_path,
            chrom    = chrom,
            pos_0    = pos_0,
            alt      = v["alt"],
            ref      = v["ref"]
        )
        if hotspot is not None:
            label  = LABEL_PCR_ARTEFACT
            reason = (
                f"Indel overlaps {hotspot.region_type} hotspot: "
                f"({hotspot.unit})x{hotspot.copies:.1f} at "
                f"{hotspot.contig}:{hotspot.start+1}-{hotspot.end} "
                f"[{hotspot.risk_level} risk]. "
                f"PCR slippage is a likely cause of this indel."
            )
        elif mean_qual >= high_quality_threshold:
            label  = LABEL_HIGH_NATURAL
            reason = (
                f"Indel outside any slippage hotspot. "
                f"Mean supporting-read Phred quality = {mean_qual:.1f} "
                f"(≥ threshold {high_quality_threshold}). "
                f"High-quality support suggests a genuine biological mutation."
            )
        else:
            label  = LABEL_LOW_NATURAL
            reason = (
                f"Indel outside any slippage hotspot. "
                f"Mean supporting-read Phred quality = {mean_qual:.1f} "
                f"(< threshold {high_quality_threshold}). "
                f"Low quality support — possible sequencing error or low-frequency artefact."
            )
        rows.append({
            "chrom":                   chrom,
            "pos_1":                   v["pos_1"],
            "ref":                     v["ref"],
            "alt":                     v["alt"],
            "mut_type":                v["mut_type"],
            "mut_size_bp":             v["mut_size"],
            "total_depth":             v["depth"],
            "ref_depth":               v["ref_depth"],
            "alt_depth":               v["alt_depth"],
            "vaf_percent":             round(v["vaf"] * 100, 2),
            "mean_alt_supporting_phred": mean_qual,
            "n_alt_supporting_reads":  n_alt_reads,
            "in_hotspot":              hotspot is not None,
            "hotspot_type":            hotspot.region_type  if hotspot else "none",
            "hotspot_unit":            hotspot.unit         if hotspot else "none",
            "hotspot_copies":          hotspot.copies       if hotspot else 0,
            "hotspot_risk_level":      hotspot.risk_level   if hotspot else "none",
            "hotspot_start_1based":    hotspot.start + 1    if hotspot else None,
            "hotspot_end_1based":      hotspot.end          if hotspot else None,
            "classification":          label,
            "classification_reason":   reason,
        })
        if logger and (i + 1) % 100 == 0:
            logger.debug(f"  Classified {i+1}/{len(variants)} variants...")
    df = pd.DataFrame(rows)
    if df.empty:
        if logger:
            logger.warning("  No variants passed all filters — output DataFrame is empty.")
        return df
    if logger:
        counts = df["classification"].value_counts()
        logger.info("  Classification summary:")
        for label, count in counts.items():
            logger.info(f"    {label:<40} : {count:,}")
    return df
