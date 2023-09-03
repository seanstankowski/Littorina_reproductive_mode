#!/usr/bin/env python

import argparse
import sys
import gzip
import ete3
import operator
import twisst
import itertools
import numpy as np

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread
from time import sleep

##############################################################################################################################

'''A function that reads from the line queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def weightTree_wrapper(lineQueue, resultQueue, taxa, taxonNames, outgroup, nIts=None,
                       topos = None, getLengths=False, getDists = False, method = "fixed", thresholdDict=None):
    nTaxa=len(taxa)
    nTopos = len(topos)
    while True:
        lineNumber,line = lineQueue.get()
        tree = twisst.readTree(line)
        
        if tree:
            #clean off unneccesary branches - speeds up downstream analyses
            leafNamesSet = set([leaf.name for leaf in tree.get_leaves()])
            if namesSet != leafNamesSet:
                assert namesSet.issubset(leafNamesSet), "Named samples not present in tree."
                tree = twisst.getPrunedCopy(tree, leavesToKeep=names, preserve_branch_length=True)
            
            if outgroup is not None: tree.set_outgroup(outgroup)
            
            if verbose: print >> sys.stderr, "Getting weights using method:", method

            weightsData = None
            
            if method == "complete":
                weightsData = twisst.weightTreeSimp(tree=tree, taxa=taxa, taxonNames=taxonNames, topos=topos, getLengths=getLengths, getDists=getDists, abortCutoff=args.abortCutoff)
            
            if method == "fixed" or (method == "complete" and backupMethod == "fixed" and weightsData == None):
                weightsData = twisst.weightTree(tree=tree, taxa=taxa, taxonNames=taxonNames, nIts=nIts, topos=topos, getLengths=getLengths, getDists=getDists)
            
            if method == "threshold" or (method == "complete" and backupMethod == "threshold" and weightsData == None):
                weightsData = twisst.weightTreeThreshold(tree=tree, taxa=taxa, taxonNames=taxonNames, thresholdDict=thresholdDict, topos=topos, getLengths=getLengths, getDists=getDists)

            weightsLine = "\t".join([str(x) for x in weightsData["weights"]])

            if getDists:
                distsByTopo = []
                for x in range(nTopos):
                    distsByTopo.append("\t".join([str(round(weightsData["dists"][pair[0],pair[1],x], 4)) for pair in itertools.combinations(range(nTaxa), 2)]))
                distsLine = "\t".join(distsByTopo)
            if getLengths:
                lengthsLine = "\t".join(["\t".join([str(round(l,4)) for l in weightsData["lengths"][x]]) for x in range(nTopos)])
        else:
            weightsLine="\t".join(["nan"]*nTopos)
            if getDists: distsLine = "\t".join(["nan"]*nTopos*len(list(itertools.combinations(range(nTaxa), 2))))
            if getLengths: lengthsLine = "\t".join(["nan"]*sum(map(len,children)))
        
        if verbose: print >> sys.stderr, "Analysed tree", lineNumber
        result = [weightsLine]
        if getDists: result.append(distsLine)
        if getLengths: result.append(lengthsLine)
        resultQueue.put((lineNumber, tuple(result), True))
        

'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(resultQueue, writeQueue, verbose):
  global resultsReceived
  sortBuffer = {}
  expect = 0
  while True:
    lineNumber,result,good = resultQueue.get()
    resultsReceived += 1
    if verbose: print >> sys.stderr, "Sorter received result", lineNumber
    if lineNumber == expect:
      writeQueue.put((lineNumber,result,good))
      if verbose: print >> sys.stderr, "Result", lineNumber, "sent to writer"
      expect +=1
      #now check buffer for further results
      while True:
        try:
          result,good = sortBuffer.pop(str(expect))
          writeQueue.put((expect,result,good))
          if verbose: print >> sys.stderr, "Result", expect, "sent to writer"
          expect +=1
        except:
          break
    else:
      #otherwise this line is ahead of us, so add to buffer dictionary
      sortBuffer[str(lineNumber)] = (result,good)

'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, outs):
    global resultsWritten
    global resultsHandled
    while True:
        lineNumber,result,good = writeQueue.get()
        if verbose:
            print >> sys.stderr, "Writer received result", lineNumber
            if good: print >> sys.stderr, "Writing good result."
            else: print >> sys.stderr, "Omitting bad result."
        if good:
            for x in range(len(outs)):
                outs[x].write(result[x] + "\n")
            resultsWritten += 1
        resultsHandled += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        print >> sys.stderr, linesQueued, "trees queued |", resultsReceived, "results received |", resultsWritten, "results written."

############################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", action = "store")
parser.add_argument("-w", "--weightsFile", help="Output file of all weights", action = "store")
parser.add_argument("-D", "--distsFile", help="Output file of mean pairwise dists", action = "store", required = False)
parser.add_argument("-L", "--lengthsFile", help="Output file of average branch lengths", action = "store", required = False)
parser.add_argument("--inputTopos", help="Input file for user-defined topologies (optional)", action = "store", required = False)
parser.add_argument("--outputTopos", help="Output file for topologies used", action = "store", required = False)
parser.add_argument("--outgroup", help="Outgroup for rooting - only affects speed", action = "store")
parser.add_argument("--method", help="Tree sampling method", choices=["fixed", "threshold", "complete"], action = "store", default = "fixed")
parser.add_argument("--backupMethod", help="Backup method if aborting complete", choices=["fixed", "threshold"], action = "store", default = "fixed")
parser.add_argument("--iterations", help="Number of iterations for fixed partial sampling", type=int, action = "store", default = 400)
parser.add_argument("--abortCutoff", help="# tips in simplified tree to abort 'complete' weighting", type=int, action = "store", default = 100000)
parser.add_argument("--thresholdTable", help="Lookup_table_for_sampling_thresholds", action = "store")
parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs="+", required = True, metavar=("name","[inds]"))
parser.add_argument("--groupsFile", help="Optional file of sample names and groups", action = "store", required = False)
parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", action = "store", type=int, default=1, required = False)
parser.add_argument("--verbose", help="Verbose output", action="store_true")


args = parser.parse_args()
#args = parser.parse_args("-n 5 -t test.trees -o test.topos.txt -w test.weights.B.csv -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

getDists  = args.distsFile is not None
getLengths = args.lengthsFile is not None

method = args.method
backupMethod = args.backupMethod
threads = args.threads

verbose = args.verbose

#################################################################################################################################
#parse taxa
assert len(args.group) >= 4, "Please specify at least four groups."
taxonNames = []
taxa = []
for g in args.group:
    taxonNames.append(g[0])
    if len(g) > 1: taxa.append(g[1].split(","))
    else: taxa.append([])

if args.groupsFile:
    with open(args.groupsFile, "r") as gf: groupDict = dict([ln.split() for ln in gf.readlines()])
    for sample in groupDict.keys():
        try: taxa[taxonNames.index(groupDict[sample])].append(sample)
        except: pass

assert min([len(t) for t in taxa]) >= 1, "Please specify at least one sample name per group."


names = [t for taxon in taxa for t in taxon]
namesSet = set(names)
assert len(names) == len(namesSet), "Each sample should only be in one group."

#get all topologies
if args.inputTopos:
    with open(args.inputTopos, "r") as tf: topos = [ete3.Tree(ln) for ln in tf.readlines()]
else: topos = twisst.allTopos(taxonNames, [])

for topo in topos: topo.set_outgroup(taxonNames[-1])

sys.stderr.write(twisst.asciiTrees(topos,5) + "\n")

nTopos = len(topos)

if args.outputTopos:
    with open(args.outputTopos, "w") as tf:
        tf.write("\n".join([t.write(format = 9) for t in topos]) + "\n")

#################################################################################################################################

# check method

if method == "fixed" or (method == "complete" and backupMethod == "fixed"):
    nIts = args.iterations
    if nIts >= np.prod([len(t) for t in taxa]):
        print >> sys.stderr, "Warning: number of iterations is equal or greater than possible combinations.\n"
        nIts = np.prod([len(t) for t in taxa])
        print >> sys.stderr, "This could be very slow. Use method 'complete' for fast(er) exhaustive sampling."
    thresholdDict = None
elif method == "threshold" or (method == "complete" and backupMethod == "threshold"):
    nIts = None
    assert args.thresholdTable, "A threshold table must be provided using argument --thresholdTable."
    thresholdTableFileName = args.thresholdTable
    with open(thresholdTableFileName) as ttf:
        thresholdDict = dict([(int(tries),int(threshold)) for line in ttf.readlines() for tries,threshold in (line.split(),)])
else:
    nIts = None
    thresholdDict = None
#################################################################################################################################
### file for weights

if args.weightsFile:
    weightsFile = gzip.open(args.weightsFile, "w") if args.weightsFile.endswith(".gz") else open(args.weightsFile, "w")
else: weightsFile = sys.stdout

for x in range(nTopos): weightsFile.write("#topo" + str(x+1) + " " + topos[x].write(format = 9) + "\n") 

weightsFile.write("\t".join(["topo" + str(x+1) for x in range(nTopos)]) + "\n")

outs = [weightsFile]


### file for lengths

if getDists:
    distsFile = gzip.open(args.distsFile, "w") if args.distsFile[-3:] == ".gz" else open(args.distsFile, "w")
    for x in range(nTopos):
        distsFile.write("\t".join(["topo" + str(x+1) + "_" + "_".join(pair) for pair in itertools.combinations(taxonNames,2)]) + "\t")
    distsFile.write("\n")
    outs += [distsFile]

if getLengths:
    lengthsFile = gzip.open(args.lengthsFile, "w") if args.lengthsFile[-3:] == ".gz" else open(args.lengthsFile, "w")
    #name internal nodes
    for t in topos: twisst.addNodeNames(t)
    parentsAndChildren = [[(n.up.name,n.name,) for n in t.traverse() if n.up is not None] for t in topos]
    children = [[pc[1] for pc in parentsAndChildren[x]] for x in range(nTopos)]
    for x in range(nTopos): lengthsFile.write("#" + "topo" + str(x+1) + "\t" + " ".join(["--".join(pair) for pair in parentsAndChildren[x]]) + "\n")
    #lengthsFile.write("\t".join(["\t".join(["topo" + str(x+1) + "_" + nodeName for nodeName in children[x]]) for x in range(nTopos)]) + "\n")
    outs += [lengthsFile]

############################################################################################################################################

#counting stat that will let keep track of how far we are
linesQueued = 0
resultsReceived = 0
resultsWritten = 0
resultsHandled = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
lineQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
for x in range(threads):
    worker = Process(target=weightTree_wrapper, args = (lineQueue, resultQueue, taxa, taxonNames, args.outgroup,
                                                        nIts, topos, getLengths, getDists, method,thresholdDict,))
    worker.daemon = True
    worker.start()
    print >> sys.stderr, "started worker", x

'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, outs,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()



############################################################################################################################################

#open tree file

if args.treeFile:
    treeFile = gzip.open(args.treeFile, "r") if args.treeFile.endswith(".gz") else open(args.treeFile, "r")
else: treeFile = sys.stdin

line = treeFile.readline()


##########################################################################################################################################

while len(line) >= 1:
    lineQueue.put((linesQueued,line.rstrip()))
    linesQueued += 1
    line = treeFile.readline()

############################################################################################################################################


### wait for queues to empty


print >> sys.stderr, "\nFinished reading trees...\n"
while resultsHandled < linesQueued:
    sleep(1)

sleep(5)

treeFile.close()
weightsFile.close()
if getDists: distsFile.close()
if getLengths: lengthsFile.close()


print >> sys.stderr, str(linesQueued), "lines were read.\n"
print >> sys.stderr, str(resultsWritten), "results were written.\n"

sys.exit()





