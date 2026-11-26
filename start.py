"""Start the same review desk straight from its unpacked folder."""
from pathlib import Path
import sys
if sys.version_info < (3,10):raise SystemExit('请使用 Python 3.10 或更新版本，例如在 Miniforge 环境中运行。')
sys.path.insert(0,str(Path(__file__).resolve().parent/'source'))
from brca_review.desk import main
if __name__=='__main__':main()
