import subprocess
import sys

from dune.testtools.parametertree.parser import parse_ini_file

if sys.argv[1] == "example":
    subprocess.call("cp /inst/*.ini .", shell=True)
    print("Copied some example configuration files into the working directory!")
    sys.exit(0)

ini = parse_ini_file(sys.argv[1])
app = ini.get("app", "cubeapp")

ret = subprocess.call(["/inst/{}".format(app), sys.argv[1]])
sys.exit(ret)
