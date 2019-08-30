import sys
from dune.testtools.parametertree.parser import parse_ini_file

ini = parse_ini_file(sys.argv[1])
sys.stdout.write(ini.get("grid.filename"))
sys.exit(0)
