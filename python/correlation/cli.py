"""CLI alias forwarding to matcorr.cli."""

from matcorr.cli import main, parse_args

__all__ = ["main", "parse_args"]

if __name__ == "__main__":
    import sys
    sys.exit(main())
