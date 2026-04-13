import sys
import os

phlodist_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "src")
)
sys.path.insert(0, phlodist_dir)

if __name__ == "__main__":
    from phylodist.main import main
    main()