#!/bin/python
import sys

def readDomMap(mapfile):
    domaindict = {}
    with open(mapfile,'r') as IN:
        for line in IN:
            linestr = line.rstrip('\n')
            ID = linestr.split(",")[0]
            if ID not in domaindict:
            	domaindict[ID] = list()
            loclist = [linestr.split(",")[1]]
            for loc in loclist:
            	for subloc in loc.split("="):
            		domaindict[ID]+= [subloc]
    IN.close()
    return domaindict

def filterFasta(domaindict):
    start = 0
    end = 0
    positions = list()
    for line in sys.stdin:
        linestr = line.rstrip('\n')
        if linestr.startswith('>'):
            header = linestr
            positions = list()
            ID = linestr.lstrip('>').split(" ")[0]
            if ID in domaindict.keys():
            	positions += domaindict[ID]
        else:
            if len(positions) > 0:
            	sequence = ""
            	for cut in positions:
            		print(">" + ID + "__" + str(cut))
            		try:
            			start = int(cut.split("-")[0]) - 1
            			end = int(cut.split("-")[1]) - 1
            		except ValueError:
            			start  = 0
            			end = int(cut.split("-")[0]) - 1
            		#sequence = sequence + linestr[start:end]
            		sequence = linestr[start:end]
            		#print(str(start) + " " + str(end) + " " + str(end - start) + str(len(linestr[start:end])))
            		print(sequence)           		

if __name__ == '__main__':
    mapfile = sys.argv[1]
    domaindict = readDomMap(mapfile)
    filterFasta(domaindict)
