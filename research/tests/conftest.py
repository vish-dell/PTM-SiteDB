from pathlib import Path
import sys

# Research scripts use imports like `from scripts...`; make the research root explicit
# for test runs launched from either the project root or the research subdirectory.
RESEARCH_ROOT = Path(__file__).resolve().parents[1]
if str(RESEARCH_ROOT) not in sys.path:
    sys.path.insert(0, str(RESEARCH_ROOT))
