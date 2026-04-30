import argparse
import sys
from pathlib import Path
from dependency_check import verify_dependencies
from aligner import align_reads
from logger import setup_logger
from hotspot_finder import build_hotspot_map
from classifier import classify_variants
from reporter import generate_all_reports
from effect_predictor import predict_effects
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass
def run_pipeline(args):
    logger = setup_logger(verbose=args.verbose)
    logger.info("Starting PCR Slippage Artefact Detection Pipeline")
    verify_dependencies(aligner=args.aligner, logger=logger, require_bcl=False)
    import datetime
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_out = Path(args.output)
    out_dir = base_out.parent / f"{base_out.name}_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    intermediates_dir = out_dir / "intermediates"
    if args.keep_intermediates:
        intermediates_dir.mkdir(exist_ok=True)
    else:
        intermediates_dir = out_dir 
    logger.info("Step 1: Aligning reads...")
    bam_path = align_reads(
        fastq1=args.fastq,
        fastq2=None,
        reference=args.reference,
        output_dir=str(intermediates_dir),
        aligner=args.aligner,
        threads=args.threads,
        logger=logger
    )
    logger.info("Step 2: Calling indels...")
    import subprocess
    vcf_path = str(intermediates_dir / "calls.vcf")
    mpileup_cmd = ["bcftools", "mpileup", "-a", "FORMAT/AD", "-Ou", "-f", args.reference, bam_path]
    call_cmd = ["bcftools", "call", "-mv", "-Ov", "-o", vcf_path]
    logger.debug(f"Running: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")
    p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    out, err = p2.communicate()
    if p2.returncode != 0:
        logger.error(f"bcftools failed: {err.decode()}")
        sys.exit(1)
    logger.info("Step 3: Building slippage hotspot map...")
    hotspot_map = build_hotspot_map(
        reference_fasta=args.reference,
        min_hp_len=args.min_hp_len,
        min_str_unit=args.min_str_unit,
        max_str_unit=args.max_str_unit,
        min_str_copies=args.min_str_copies,
        output_dir=str(intermediates_dir),
        logger=logger
    )
    logger.info("Step 4: Classifying variants...")
    classified_df = classify_variants(
        vcf_path=vcf_path,
        bam_path=bam_path,
        hotspot_map=hotspot_map,
        high_quality_threshold=args.high_quality_threshold,
        min_quality=args.min_quality,
        min_vaf=args.min_vaf,
        logger=logger
    )
    classified_df = predict_effects(
        df=classified_df,
        reference_fasta=args.reference,
        window_size=500,
        logger=logger
    )
    logger.info("Step 5: Generating reports...")
    reports_dir = out_dir / "reports"
    generate_all_reports(
        classified_df=classified_df,
        hotspot_map=hotspot_map,
        reference=args.reference,
        output_dir=str(reports_dir),
        logger=logger
    )
    logger.info("Pipeline completed successfully!")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCR Slippage Artefact Detection Pipeline")
    parser.add_argument("--fastq", required=True, help="Single-end FASTQ file")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--min-hp-len", type=int, default=6)
    parser.add_argument("--min-str-unit", type=int, default=2)
    parser.add_argument("--max-str-unit", type=int, default=6)
    parser.add_argument("--min-str-copies", type=int, default=3)
    parser.add_argument("--min-quality", type=int, default=20)
    parser.add_argument("--high-quality-threshold", type=int, default=30)
    parser.add_argument("--min-depth", type=int, default=5)
    parser.add_argument("--min-vaf", type=float, default=0.01)
    parser.add_argument("--aligner", default="bwa")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--keep-intermediates", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    print("\n--- Interactive Parameter Setup ---")
    print("Press Enter to keep the default value in brackets.")
    def prompt_int(name, default):
        val = input(f"{name} [{default}]: ").strip()
        if not val:
            return default
        try:
            return int(val)
        except ValueError:
            print(f"Invalid integer. Using default: {default}")
            return default
    def prompt_float(name, default):
        val = input(f"{name} [{default}]: ").strip()
        if not val:
            return default
        try:
            return float(val)
        except ValueError:
            print(f"Invalid float. Using default: {default}")
            return default
    args.min_hp_len = prompt_int("Minimum homopolymer length", args.min_hp_len)
    args.min_str_unit = prompt_int("Minimum STR unit size", args.min_str_unit)
    args.max_str_unit = prompt_int("Maximum STR unit size", args.max_str_unit)
    args.min_str_copies = prompt_int("Minimum STR copies", args.min_str_copies)
    args.min_quality = prompt_int("Minimum base quality", args.min_quality)
    args.high_quality_threshold = prompt_int("High quality threshold", args.high_quality_threshold)
    args.min_depth = prompt_int("Minimum read depth", args.min_depth)
    args.min_vaf = prompt_float("Minimum variant allele frequency (VAF)", args.min_vaf)
    print("-----------------------------------\n")
    run_pipeline(args)
