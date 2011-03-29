import sys
import os
from common import *
script=sys.argv[1]
args=' '.join(sys.argv[2:])
ossystem(format('R {args} --no-save --args --slave <{script}',locals()))
