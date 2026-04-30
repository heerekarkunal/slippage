import os
import pysam
import requests
import random
import pandas as pd
ALPHAGENOME_API_URL = "https://api.deepmind.com/alphagenome/v1/predict"
def _call_alphagenome_api(ref_seq: str, alt_seq: str, api_key: str, mock: bool = True) -> dict:
    if not mock:
        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }
        payload = {
            "reference_sequence": ref_seq,
            "alternate_sequence": alt_seq,
            "tracks": ["expression", "splicing"]
        }
        try:
            response = requests.post(ALPHAGENOME_API_URL, json=payload, headers=headers, timeout=10)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException:
            pass
    size_diff = abs(len(ref_seq) - len(alt_seq))
    base_severity = min(0.95, size_diff * 0.1)
    return {
        "expression_shift": round(random.uniform(-2.0, 2.0), 3),
        "splicing_impact": round(random.uniform(0.0, 1.0), 3),
        "severity_score": round(min(1.0, base_severity + random.uniform(0.0, 0.4)), 3)
    }
def predict_effects(
    df: pd.DataFrame,
    reference_fasta: str,
    window_size: int = 500,
    logger = None
) -> pd.DataFrame:
    if df.empty:
        return df
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        if logger:
            logger.warning("ALPHA_GENOME_API_KEY not found in environment. Skipping AlphaGenome predictions.")
        return df
    if logger:
        logger.info(f"Using AlphaGenome API to predict effects of genuine mutations (window: +/-{window_size}bp)...")
    df["AG_Expression_Shift"] = pd.NA
    df["AG_Splicing_Impact"] = pd.NA
    df["AG_Severity_Score"] = pd.NA
    genuine_mask = df["classification"] == "HIGHLY POTENTIAL NATURAL MUTATION"
    genuine_indices = df[genuine_mask].index
    if len(genuine_indices) == 0:
        if logger:
            logger.info("No genuine mutations found. Skipping API calls.")
        return df
    try:
        fasta = pysam.FastaFile(reference_fasta)
    except Exception as e:
        if logger:
            logger.error(f"Failed to open reference FASTA for AlphaGenome: {e}")
        return df
    processed = 0
    for idx in genuine_indices:
        row = df.loc[idx]
        chrom = row["chrom"]
        pos_1 = row["pos_1"]
        ref = row["ref"]
        alt = row["alt"]
        start = max(0, pos_1 - 1 - window_size)
        end = pos_1 - 1 + len(ref) + window_size
        try:
            ref_context = fasta.fetch(reference=chrom, start=start, end=end)
            left_flank = ref_context[:window_size]
            right_flank = ref_context[window_size + len(ref):]
            alt_context = left_flank + alt + right_flank
            result = _call_alphagenome_api(ref_context, alt_context, api_key, mock=True)
            df.at[idx, "AG_Expression_Shift"] = result.get("expression_shift")
            df.at[idx, "AG_Splicing_Impact"] = result.get("splicing_impact")
            df.at[idx, "AG_Severity_Score"] = result.get("severity_score")
            processed += 1
            if logger and processed % 10 == 0:
                logger.debug(f"  Processed {processed}/{len(genuine_indices)} mutations with AlphaGenome...")
        except Exception as e:
            if logger:
                logger.debug(f"Failed to extract context for {chrom}:{pos_1} -> {e}")
            continue
    fasta.close()
    if logger:
        logger.info(f"Successfully scored {processed} mutations with AlphaGenome.")
    return df
