import shutil
import subprocess
import sys

REQUIRED_PYTHON_PACKAGES = {
    "pysam":      "pip install pysam",
    "Bio":        "pip install biopython",
    "pandas":     "pip install pandas",
    "numpy":      "pip install numpy",
    "matplotlib": "pip install matplotlib",
    "tqdm":       "pip install tqdm",
    "regex":      "pip install regex",
}

REQUIRED_SYSTEM_TOOLS = {
    "samtools": "sudo apt install samtools  OR  conda install -c bioconda samtools",
    "bcftools": "sudo apt install bcftools  OR  conda install -c bioconda bcftools",
    "bgzip":    "sudo apt install tabix      OR  conda install -c bioconda tabix",
    "tabix":    "sudo apt install tabix      OR  conda install -c bioconda tabix",
}

ALIGNER_TOOLS = {
    "bwa-mem2": "conda install -c bioconda bwa-mem2",
    "bwa":      "conda install -c bioconda bwa  OR  sudo apt install bwa",
    "minimap2": "conda install -c bioconda minimap2  OR  pip install minimap2",
}


def _check_python_package(package: str) -> bool:
    try:
        __import__(package)
        return True
    except ImportError:
        return False


def _check_system_tool(tool: str) -> bool:
    return shutil.which(tool) is not None


def _get_tool_version(tool: str) -> str:
    try:
        result = subprocess.run(
            [tool, "--version"],
            capture_output=True, text=True, timeout=5,
        )
        out = (result.stdout + result.stderr).strip()
        return out.split("\n")[0][:60]
    except Exception:
        return "version unknown"


def verify_dependencies(aligner: str, logger):
    missing = []

    logger.debug("Checking Python packages...")
    for pkg, install_hint in REQUIRED_PYTHON_PACKAGES.items():
        if _check_python_package(pkg):
            logger.debug(f"  ✓ Python: {pkg}")
        else:
            logger.error(f"  ✗ Python package missing: {pkg}  →  {install_hint}")
            missing.append(pkg)

    logger.debug("Checking system tools...")
    for tool, install_hint in REQUIRED_SYSTEM_TOOLS.items():
        if _check_system_tool(tool):
            ver = _get_tool_version(tool)
            logger.debug(f"  ✓ {tool}  ({ver})")
        else:
            logger.error(f"  ✗ System tool missing: {tool}  →  {install_hint}")
            missing.append(tool)

    if aligner in ALIGNER_TOOLS:
        if _check_system_tool(aligner):
            ver = _get_tool_version(aligner)
            logger.debug(f"  ✓ {aligner}  ({ver})")
        else:
            hint = ALIGNER_TOOLS[aligner]
            logger.error(f"  ✗ Aligner missing: {aligner}  →  {hint}")
            missing.append(aligner)

    if missing:
        logger.error("")
        logger.error("The following dependencies are missing:")
        for m in missing:
            logger.error(f"   - {m}")
        logger.error("")
        logger.error("Install them and re-run. See messages above for install hints.")
        sys.exit(1)

    logger.info("All dependencies satisfied.")
