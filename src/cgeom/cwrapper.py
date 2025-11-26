from pathlib import Path
import tinycwrap

base_path = Path(__file__).parent/"clib"

clib = tinycwrap.CModule(base_path/"base.c", base_path/"path.c")