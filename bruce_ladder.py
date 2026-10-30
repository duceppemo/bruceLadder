#!/usr/bin/env python3

"""
Kept so that `python bruce_ladder.py` keeps working from a clone, which is how every example in
the README and the wiki invokes it, and what the bruce_rates.py and suis_rates.py workflows
expect. The code lives in the brucetyper package; `brucetyper` is the installed entry point.
"""

from brucetyper.cli import main

__author__ = 'duceppemo'


if __name__ == '__main__':
    main()
