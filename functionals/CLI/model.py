# External
from typing import Any, Dict as D
from os.path import join
# Internal Modules
from dbgen import parser, ConnectInfo
from functionals.model import make_model
"""
Run the model defined in /functionals/model.
"""
##############################################################################


root = '/'+join(*__file__.split('/')[:-3], 'data/')


def main(args: D[str, Any]) -> None:
    """Run the model with no extensi/ons from command line."""

    m = make_model()
    db = ConnectInfo.from_file(root+'functionals.json')
    mdb = ConnectInfo.from_file(root+'functionals_log.json')
    m.run(db, mdb, **args)


if __name__ == '__main__':
    args = parser.parse_args()
    main(vars(args))
