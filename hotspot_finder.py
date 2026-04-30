import re
import regex
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from pathlib import Path
from Bio import SeqIO
@dataclass
class HotspotRegion:
    contig:       str
    start:        int           
    end:          int           
    total_length: int
    region_type:  str           
    unit:         str           
    unit_length:  int           
    copies:       float         
    risk_level:   str           
    description:  str = field(init=False)
    def __post_init__(self):
        self.description = (
            f"{self.region_type.upper()} ({self.unit})x{self.copies:.1f} "
            f"len={self.total_length} [{self.risk_level}]"
        )
    def contains(self, position: int) -> bool:
        return self.start <= position < self.end
    def overlaps(self, pos_start: int, pos_end: int) -> bool:
        return self.start < pos_end and self.end > pos_start
def _tier_risk(total_length: int, region_type: str) -> str:
    if region_type == "homopolymer":
        return "HIGH" if total_length > 12 else "MODERATE"
    else:  
        return "HIGH" if total_length > 16 else "MODERATE"
def _find_homopolymers(contig: str, sequence: str, min_len: int) -> List[HotspotRegion]:
    hits = []
    pattern = re.compile(r'([ACGTNacgtn])\1{' + str(min_len - 1) + r',}')
    for match in pattern.finditer(str(sequence)):
        seq_str = match.group(0).upper()
        unit    = match.group(1).upper()
        length  = len(seq_str)
        hits.append(HotspotRegion(
            contig       = contig,
            start        = match.start(),
            end          = match.end(),
            total_length = length,
            region_type  = "homopolymer",
            unit         = unit,
            unit_length  = 1,
            copies       = float(length),
            risk_level   = _tier_risk(length, "homopolymer")
        ))
    return hits
def _find_strs(
    contig: str,
    sequence: str,
    min_unit: int,
    max_unit: int,
    min_copies: int
) -> List[HotspotRegion]:
    hits = []
    seq_upper = str(sequence).upper()
    for unit_len in range(min_unit, max_unit + 1):
        pattern = regex.compile(
            r'([ACGTN]{' + str(unit_len) + r'})\1{' + str(min_copies - 1) + r',}',
            flags=regex.IGNORECASE
        )
        for match in pattern.finditer(seq_upper):
            unit       = match.group(1).upper()
            full_match = match.group(0)
            length     = len(full_match)
            copies     = length / unit_len
            if len(set(unit)) == 1:
                continue
            if all(b == 'N' for b in unit):
                continue
            hits.append(HotspotRegion(
                contig       = contig,
                start        = match.start(),
                end          = match.end(),
                total_length = length,
                region_type  = "STR",
                unit         = unit,
                unit_length  = unit_len,
                copies       = copies,
                risk_level   = _tier_risk(length, "STR")
            ))
    return hits
def _merge_overlapping(regions: List[HotspotRegion]) -> List[HotspotRegion]:
    if not regions:
        return []
    regions = sorted(regions, key=lambda r: (r.start, r.end))
    merged  = [regions[0]]
    for current in regions[1:]:
        last = merged[-1]
        if current.start < last.end:                        
            if current.end > last.end:
                new_end    = current.end
                new_len    = new_end - last.start
                new_risk   = "HIGH" if "HIGH" in (last.risk_level, current.risk_level) else "MODERATE"
                merged[-1] = HotspotRegion(
                    contig       = last.contig,
                    start        = last.start,
                    end          = new_end,
                    total_length = new_len,
                    region_type  = last.region_type,
                    unit         = last.unit,
                    unit_length  = last.unit_length,
                    copies       = new_len / max(last.unit_length, 1),
                    risk_level   = new_risk
                )
        else:
            merged.append(current)
    return merged
def build_hotspot_map(
    reference_fasta: str,
    min_hp_len:    int = 6,
    min_str_unit:  int = 2,
    max_str_unit:  int = 6,
    min_str_copies:int = 3,
    output_dir:    str = ".",
    logger         = None
) -> Dict[str, List[HotspotRegion]]:
    ref_path = Path(reference_fasta)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    hotspot_map: Dict[str, List[HotspotRegion]] = {}
    total_regions = 0
    total_high    = 0
    if logger:
        logger.info(f"  Scanning reference: {ref_path.name}")
        logger.info(f"  Homopolymer min length : {min_hp_len} bp")
        logger.info(f"  STR unit range         : {min_str_unit}–{max_str_unit} bp")
        logger.info(f"  STR min copies         : {min_str_copies}")
    for record in SeqIO.parse(str(ref_path), "fasta"):
        contig   = record.id
        sequence = str(record.seq).upper()
        if logger:
            logger.info(f"  Scanning contig: {contig} ({len(sequence):,} bp)")
        hp_regions  = _find_homopolymers(contig, sequence, min_len=min_hp_len)
        str_regions = _find_strs(contig, sequence, min_str_unit, max_str_unit, min_str_copies)
        all_regions = _merge_overlapping(hp_regions + str_regions)
        hotspot_map[contig] = all_regions
        n_high = sum(1 for r in all_regions if r.risk_level == "HIGH")
        total_regions += len(all_regions)
        total_high    += n_high
        if logger:
            logger.info(
                f"    Found {len(all_regions)} hotspot region(s)  "
                f"({n_high} HIGH-risk, {len(all_regions)-n_high} MODERATE)"
            )
    bed_path = out_path / "hotspots.bed"
    _write_bed(hotspot_map, bed_path)
    tsv_path = out_path / "hotspots.tsv"
    _write_tsv(hotspot_map, tsv_path)
    if logger:
        logger.info(f"  Total hotspot regions: {total_regions}  ({total_high} HIGH-risk)")
        logger.info(f"  Hotspot BED file : {bed_path}")
        logger.info(f"  Hotspot TSV file : {tsv_path}")
    return hotspot_map
def _write_bed(hotspot_map, bed_path):
    with open(bed_path, "w") as fh:
        fh.write("track name='Slippage Hotspots' description='Homopolymers and STRs'\n")
        for contig, regions in hotspot_map.items():
            for r in regions:
                score = 1000 if r.risk_level == "HIGH" else 500
                fh.write(f"{r.contig}\t{r.start}\t{r.end}\t{r.description}\t{score}\t+\n")
def _write_tsv(hotspot_map, tsv_path):
    with open(tsv_path, "w") as fh:
        fh.write(
            "contig\tstart_1based\tend_1based\ttotal_length\t"
            "region_type\tunit\tunit_length\tcopies\trisk_level\n"
        )
        for contig, regions in hotspot_map.items():
            for r in regions:
                fh.write(
                    f"{r.contig}\t{r.start+1}\t{r.end}\t{r.total_length}\t"
                    f"{r.region_type}\t{r.unit}\t{r.unit_length}\t{r.copies:.1f}\t"
                    f"{r.risk_level}\n"
                )
def position_in_hotspot(
    contig: str,
    position: int,
    hotspot_map: Dict[str, List[HotspotRegion]]
) -> Optional[HotspotRegion]:
    regions = hotspot_map.get(contig, [])
    for region in regions:
        if region.contains(position):
            return region
    return None
