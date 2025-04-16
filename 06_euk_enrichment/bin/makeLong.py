#!/bin/python
import sys

def makeLong():
	for line in sys.stdin:
		linestr = line.rstrip('\n')
		linelist = linestr.split('\t')
		if linelist[3] != "":
			uniref50 = linelist[0]
			memberlist = linelist[3].split(",")
			for member in memberlist:
				print(uniref50 + "\t" + member)

if __name__ == '__main__':
	makeLong()
