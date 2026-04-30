import subprocess
from pathlib import Path
from typing import Optional
def align_reads(
    fastq1: str,
    fastq2: Optional[str],
    reference: str,
    output_dir: str,
    aligner: str = "bwa-mem2",
    threads: int = 4,
    logger=None
) -> str:
    out     = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ref     = Path(reference)
    sample  = Path(fastq1).name.split(".fastq")[0].replace("_R1_001", "").replace("_R1", "")
    bam_out = out / f"{sample}.sorted.bam"
    if logger:
        logger.info(f"  Aligner   : {aligner}")
        logger.info(f"  Reference : {ref.name}")
        logger.info(f"  Reads     : {Path(fastq1).name}" +
                    (f" + {Path(fastq2).name}" if fastq2 else " (single-end)"))
        logger.info(f"  Output BAM: {bam_out.name}")
    _index_reference(reference, aligner, threads, logger)
    sam_path = out / f"{sample}.sam"
    if aligner in ("bwa-mem2", "bwa"):
        _run_bwa(aligner, reference, fastq1, fastq2, sam_path, threads, sample, logger)
    elif aligner == "minimap2":
        _run_minimap2(reference, fastq1, fastq2, sam_path, threads, logger)
    _sam_to_sorted_bam(sam_path, bam_out, threads, logger)
    sam_path.unlink(missing_ok=True)  
    _run(["samtools", "index", str(bam_out)], "Indexing BAM", logger)
    stats_path = out / f"{sample}.flagstat.txt"
    with open(stats_path, "w") as fh:
        subprocess.run(["samtools", "flagstat", str(bam_out)], stdout=fh)
    if logger:
        _log_flagstat(stats_path, logger)
    return str(bam_out)
def _index_reference(reference: str, aligner: str, threads: int, logger):
    ref = Path(reference)
    if aligner == "bwa-mem2":
        index_file = Path(str(reference) + ".bwt.2bit.64")
        if index_file.exists():
            if logger: logger.debug("  Reference index (bwa-mem2) already exists — skipping")
            return
        if logger: logger.info("  Indexing reference with bwa-mem2...")
        _run(["bwa-mem2", "index", str(ref)], "bwa-mem2 index", logger)
    elif aligner == "bwa":
        index_file = Path(str(reference) + ".bwt")
        if index_file.exists():
            if logger: logger.debug("  Reference index (bwa) already exists — skipping")
            return
        if logger: logger.info("  Indexing reference with bwa...")
        _run(["bwa", "index", str(ref)], "bwa index", logger)
    elif aligner == "minimap2":
        mmi = Path(str(reference) + ".mmi")
        if mmi.exists():
            if logger: logger.debug("  Reference index (minimap2) already exists — skipping")
            return
        if logger: logger.info("  Indexing reference with minimap2...")
        _run(["minimap2", "-d", str(mmi), str(ref)], "minimap2 index", logger)
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        _run(["samtools", "faidx", str(ref)], "samtools faidx", logger)
def _run_bwa(aligner, reference, fastq1, fastq2, sam_out, threads, sample, logger):
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:lib1"
    cmd = [
        aligner, "mem",
        "-t", str(threads),
        "-R", rg,
        str(reference),
        str(fastq1),
    ]
    if fastq2:
        cmd.append(str(fastq2))
    if logger:
        logger.info("  Running alignment (this may take a few minutes)...")
        logger.debug(f"  CMD: {' '.join(cmd[:6])} ...")
    with open(sam_out, "w") as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{aligner} alignment failed:\n{result.stderr}")
def _run_minimap2(reference, fastq1, fastq2, sam_out, threads, logger):
    cmd = [
        "minimap2",
        "-ax", "sr",          
        "-t", str(threads),
        "--secondary=no",
        str(reference),
        str(fastq1),
    ]
    if fastq2:
        cmd.append(str(fastq2))
    if logger:
        logger.info("  Running minimap2 alignment...")
        logger.debug(f"  CMD: {' '.join(cmd[:6])} ...")
    with open(sam_out, "w") as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 alignment failed:\n{result.stderr}")
def _sam_to_sorted_bam(sam_path: Path, bam_out: Path, threads: int, logger):
    if logger:
        logger.info("  Converting SAM → sorted BAM...")
    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(bam_out),
        str(sam_path)
    ]
    result = subprocess.run(sort_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools sort failed:\n{result.stderr}")
def _log_flagstat(stats_path: Path, logger):
    try:
        lines = stats_path.read_text().strip().split("\n")
        for line in lines[:6]:
            logger.info(f"  {line.strip()}")
    except Exception:
        pass
def _run(cmd: list, label: str, logger):
    if logger:
        logger.debug(f"  Running: {label}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr}")
