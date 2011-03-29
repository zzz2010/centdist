import sys
import os
from common import *
script=sys.argv[1]
args=' '.join(map(lambda x:"'"+x+"'",sys.argv[2:]))
ossystem(format('R --slave --no-save --args {args} <{script}',locals()))
