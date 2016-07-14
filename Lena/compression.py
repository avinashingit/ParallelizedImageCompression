import time
import heapq
import os
import math
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
from PIL import Image
from multiprocessing import Process
import numpy as np
np.set_printoptions(threshold='nan')

def nob(n):
    return int(math.log(n, 2)) + 1

def doubling_range(start, stop):
    while start < stop:
        yield start
        start <<= 1

class Node(object):
    def __init__(self, pairs, frequency):
        # print pairs
        self.pairs = pairs
        self.frequency = frequency

    def __repr__(self):
        return repr(self.pairs) + ", " + repr(self.frequency)

    def merge(self, other):
        total_frequency = self.frequency + other.frequency
        for p in self.pairs:
            p[1] = "1" + p[1]
        for p in other.pairs:
            p[1] = "0" + p[1]
        new_pairs = self.pairs + other.pairs
        return Node(new_pairs, total_frequency)

    def __lt__(self, other):
        return self.frequency < other.frequency

def imread(imageFile):
        #read image
        currentImage = Image.open(imageFile)
        image = currentImage.load()
        [rows,columns] = currentImage.size
        redFile = open("redFile.txt","w")
        greenFile = open("greenFile.txt","w")
        blueFile = open("blueFile.txt","w")
        for i in range(0,rows):
            redFile.write("\n")
            greenFile.write("\n")
            blueFile.write("\n")
            for j in range(0, columns):
                redFile.write(str(image[i,j][0])+" ")
                greenFile.write(str(image[i,j][1])+" ")
                blueFile.write(str(image[i,j][2])+" ")

def getImageSize(imageFile):
        currentImage = Image.open(imageFile)
        return currentImage.size

def blockshaped(arr, nrows, ncols):
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def unblockshaped(arr, h, w):
    n, nrows, ncols = arr.shape
    return (arr.reshape(h//nrows, -1, nrows, ncols)
               .swapaxes(1,2)
               .reshape(h, w))

def clustering(i,string,y):
    z = y[i].reshape(1,blockSize*blockSize)
    centroids,_ = kmeans(z[0],numberOfClusters)
    # print centroids
    k = 0;
    redCluster = open(string+"Cluster"+str(i)+".txt","w")
    with open(string+"ClusterTable"+str(i)+".txt","w") as redClusterTable:
        for identifier in centroids:
            redClusterTable.write(str(k))
            redClusterTable.write(" "+str(identifier)+"\n");
            k = k+1
    idx,_ = vq(z[0],centroids)
    c = 1
    a = list()
    for cid in idx:
        redCluster.write(str(cid)+" ")
        a.append(cid)
        c = c+1
        if(c==blockSize+1):
            redCluster.write("\n")
            c = 1;

def runParallelClustering(string):
    red = []
    with open(string+"File.txt") as redFile:
        y = redFile.read()
        for m in y.split():
            red.append(int(m))
    redarray  = np.asarray(red);
    redarray = redarray.reshape(rows,columns)
    y = blockshaped(redarray, blockSize, blockSize)
    proc = []
    for i in range(0,numberOfBlocks):
        p = Process(target=clustering, args=(i,string,y,))
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

def runAllClusteringParallel():
    proc = []
    for string in ["red","green","blue"]:
        p = Process(target=runParallelClustering, args=(string,))
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

def clusterEncoding(string):
    final = []
    for i in range(0,numberOfBlocks):
        with open(string+"Cluster"+str(i)+".txt","r") as s:
            i = s.read()
            k = []
            for m in i.split():
                k.append(int(m))
        k = np.asarray(k)
        k = k.reshape(blockSize,blockSize)
        final.append(k)
    final = np.asarray(final)
    p = unblockshaped(final,rows,columns)
    with open(string+"ClusterEncoded.txt","w") as w:
        for i in range(0,rows):
            for m in p[i]:
                w.write(str(m)+" ")
            w.write("\n")

def runClusterEncodingParallel():
    proc = []
    for string in ["red","green","blue"]:
        p = Process(target=clusterEncoding, args=(string,))
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

def miner():
    redFile = []
    with open("redClusterEncoded.txt", "r") as r:
        for line in r:
            line = line.replace("\n","")
            redFile.append(" "+line+" ")
    # print redFile
    #Find all one length sequences
    frequentPatterns = {}
    currentLength = 1
    frequentPatterns[currentLength] = {}
    for i in range(0,numberOfClusters):
        cKey = " "+str(i)+" "
        for eachLine in redFile:
            if(cKey in eachLine):
                if(cKey in frequentPatterns[currentLength].keys()):
                    frequentPatterns[currentLength][cKey] = frequentPatterns[currentLength][cKey] + 1
                else:
                    frequentPatterns[currentLength][cKey] = 1
    # print frequentPatterns
    #Enumerate length 2 frequent sequences
    currentLength = 2
    frequentPatterns[currentLength] = {}
    for key1 in frequentPatterns[currentLength-1].keys():
        for key2 in frequentPatterns[currentLength-1].keys():
            newKey = key1.rstrip()+key2
            frequentPatterns[currentLength][newKey] = 0
    for eachKey in frequentPatterns[currentLength].keys():
        for line in redFile:
            if(eachKey in line):
                frequentPatterns[currentLength][eachKey] = frequentPatterns[currentLength][eachKey] + 1
    for key in frequentPatterns[currentLength].keys():
        if(frequentPatterns[currentLength][key]<minimumSupport):
            frequentPatterns[currentLength].pop(key, None)
    currentLength = currentLength + 1
    # print frequentPatterns
    #Enumerate all frequent sequences
    for i in range(currentLength, columns+1):
        # print i
        # if(len(frequentPatterns[i].keys())!=0):
        frequentPatterns[i] = {}
        #Generate the i-length keys
        for key1 in frequentPatterns[i-1].keys():
            for key2 in frequentPatterns[i-1].keys():
                mkey1 = key1.lstrip().rstrip().split(" ")
                mkey2 = key2.lstrip().rstrip().split(" ")
                if(mkey1[1:len(mkey1)]==mkey2[0:len(mkey2)-1]):
                    key2Split = key2.lstrip().rstrip().split(" ")
                    newKey = key1+key2Split[len(key2Split)-1]+" "
                    frequentPatterns[i][newKey] = 0;
        # print "Loop"
        #Find the support of each i-length key
        for line in redFile:
            for eachKey in frequentPatterns[i].keys():
                if(eachKey in line):
                    frequentPatterns[i][eachKey] = frequentPatterns[i][eachKey] + 1

        #Pop all sequences that don't satisfy the minimum support threshold
        for key in frequentPatterns[i].keys():
            if(frequentPatterns[i][key]<minimumSupport):
                frequentPatterns[i].pop(key, None)

        #Pop all non closed sequences
        for key1 in frequentPatterns[i].keys():
            for key2 in frequentPatterns[i-1].keys():
                if((key2 in key1) and (frequentPatterns[i][key1] == frequentPatterns[i-1][key2])):
                    frequentPatterns[i-1].pop(key2, None)

    allKeys = []
    allFrequentPatterns = {}
    for i in range(1, columns+1):
        for key in frequentPatterns[i].keys():
            allKeys.insert(0,key.lstrip().rstrip())
            allFrequentPatterns[key.lstrip().rstrip()] = 0
            # print key+"---"+str(frequentPatterns[i][key])

    #By this time all frequent sequences are generated
    '''*********************************************'''
    #Find the modified frequency of each sequence
    for line in redFile:
        # print "x"
        for key in allKeys:
            # print key
            if(len(key.split(" ")) == 1):
                # print key
                key = " "+key+" "
                if(key in line):
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
                    # print line
            else:
                key = " "+key+" "
                if(key in line):
                    # key = key.lstrip().rstrip()
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
        # print line
    finalKeys = open("redPatterns.txt","w")
    finalLength = 0
    for eachKey in allKeys:
        if(allFrequentPatterns[eachKey]!=0):
            finalKeys.write(eachKey.lstrip().rstrip()+"-"+str(allFrequentPatterns[eachKey])+"\n")
            finalLength = finalLength + (len(eachKey.lstrip().rstrip().split(" "))*allFrequentPatterns[eachKey])


    greenFile = []
    with open("greenClusterEncoded.txt", "r") as r:
        for line in r:
            line = line.replace("\n","")
            greenFile.append(" "+line+" ")

    #Find all one length sequences
    frequentPatterns = {}
    currentLength = 1
    frequentPatterns[currentLength] = {}
    for i in range(0,numberOfClusters):
        cKey = " "+str(i)+" "
        for eachLine in greenFile:
            if(cKey in eachLine):
                if(cKey in frequentPatterns[currentLength].keys()):
                    frequentPatterns[currentLength][cKey] = frequentPatterns[currentLength][cKey] + 1
                else:
                    frequentPatterns[currentLength][cKey] = 1

    #Enumerate length 2 frequent sequences
    currentLength = 2
    frequentPatterns[currentLength] = {}
    for key1 in frequentPatterns[currentLength-1].keys():
        for key2 in frequentPatterns[currentLength-1].keys():
            newKey = key1.rstrip()+key2
            frequentPatterns[currentLength][newKey] = 0
    for eachKey in frequentPatterns[currentLength].keys():
        for line in redFile:
            if(eachKey in line):
                frequentPatterns[currentLength][eachKey] = frequentPatterns[currentLength][eachKey] + 1
    for key in frequentPatterns[currentLength].keys():
        if(frequentPatterns[currentLength][key]<minimumSupport):
            frequentPatterns[currentLength].pop(key, None)
    currentLength = currentLength + 1

    #Enumerate all frequent sequences
    for i in range(currentLength, columns+1):
        # if(len(frequentPatterns[i].keys())!=0):
        frequentPatterns[i] = {}
        #Generate the i-length keys
        for key1 in frequentPatterns[i-1].keys():
            for key2 in frequentPatterns[i-1].keys():
                mkey1 = key1.lstrip().rstrip().split(" ")
                mkey2 = key2.lstrip().rstrip().split(" ")
                if(mkey1[1:len(mkey1)]==mkey2[0:len(mkey2)-1]):
                    key2Split = key2.lstrip().rstrip().split(" ")
                    newKey = key1+key2Split[len(key2Split)-1]+" "
                    frequentPatterns[i][newKey] = 0;

        #Find the support of each i-length key
        for line in greenFile:
            for eachKey in frequentPatterns[i].keys():
                if(eachKey in line):
                    frequentPatterns[i][eachKey] = frequentPatterns[i][eachKey] + 1

        #Pop all sequences that don't satisfy the minimum support threshold
        for key in frequentPatterns[i].keys():
            if(frequentPatterns[i][key]<minimumSupport):
                frequentPatterns[i].pop(key, None)

        #Pop all non closed sequences
        for key1 in frequentPatterns[i].keys():
            for key2 in frequentPatterns[i-1].keys():
                if((key2 in key1) and (frequentPatterns[i][key1] == frequentPatterns[i-1][key2])):
                    frequentPatterns[i-1].pop(key2, None)

    allKeys = []
    allFrequentPatterns = {}
    for i in range(1, columns+1):
        for key in frequentPatterns[i].keys():
            allKeys.insert(0,key.lstrip().rstrip())
            allFrequentPatterns[key.lstrip().rstrip()] = 0
            # print key+"---"+str(frequentPatterns[i][key])

    #By this time all frequent sequences are generated
    '''*********************************************'''
    #Find the modified frequency of each sequence
    for line in greenFile:
        for key in allKeys:
            if(len(key.split(" ")) == 1):
                # print key
                key = " "+key+" "
                if(key in line):
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
                    # print line
            else:
                key = " "+key+" "
                if(key in line):
                    # key = key.lstrip().rstrip()
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
    finalKeys = open("greenPatterns.txt","w")
    finalLength = 0
    for eachKey in allKeys:
        if(allFrequentPatterns[eachKey]!=0):
            finalKeys.write(eachKey.lstrip().rstrip()+"-"+str(allFrequentPatterns[eachKey])+"\n")
            finalLength = finalLength + (len(eachKey.lstrip().rstrip().split(" "))*allFrequentPatterns[eachKey])

    blueFile = []
    with open("blueClusterEncoded.txt", "r") as r:
        for line in r:
            line = line.replace("\n","")
            blueFile.append(" "+line+" ")

    #Find all one length sequences
    frequentPatterns = {}
    currentLength = 1
    frequentPatterns[currentLength] = {}
    for i in range(0,numberOfClusters):
        cKey = " "+str(i)+" "
        for eachLine in blueFile:
            if(cKey in eachLine):
                if(cKey in frequentPatterns[currentLength].keys()):
                    frequentPatterns[currentLength][cKey] = frequentPatterns[currentLength][cKey] + 1
                else:
                    frequentPatterns[currentLength][cKey] = 1

    #Enumerate length 2 frequent sequences
    currentLength = 2
    frequentPatterns[currentLength] = {}
    for key1 in frequentPatterns[currentLength-1].keys():
        for key2 in frequentPatterns[currentLength-1].keys():
            newKey = key1.rstrip()+key2
            frequentPatterns[currentLength][newKey] = 0
    for eachKey in frequentPatterns[currentLength].keys():
        for line in redFile:
            if(eachKey in line):
                frequentPatterns[currentLength][eachKey] = frequentPatterns[currentLength][eachKey] + 1
    for key in frequentPatterns[currentLength].keys():
        if(frequentPatterns[currentLength][key]<minimumSupport):
            frequentPatterns[currentLength].pop(key, None)
    currentLength = currentLength + 1

    #Enumerate all frequent sequences
    for i in range(currentLength, columns+1):
        # if(len(frequentPatterns[i].keys())!=0):
        frequentPatterns[i] = {}
        #Generate the i-length keys
        for key1 in frequentPatterns[i-1].keys():
            for key2 in frequentPatterns[i-1].keys():
                mkey1 = key1.lstrip().rstrip().split(" ")
                mkey2 = key2.lstrip().rstrip().split(" ")
                if(mkey1[1:len(mkey1)]==mkey2[0:len(mkey2)-1]):
                    key2Split = key2.lstrip().rstrip().split(" ")
                    newKey = key1+key2Split[len(key2Split)-1]+" "
                    frequentPatterns[i][newKey] = 0;

        #Find the support of each i-length key
        for line in blueFile:
            for eachKey in frequentPatterns[i].keys():
                if(eachKey in line):
                    frequentPatterns[i][eachKey] = frequentPatterns[i][eachKey] + 1

        #Pop all sequences that don't satisfy the minimum support threshold
        for key in frequentPatterns[i].keys():
            if(frequentPatterns[i][key]<minimumSupport):
                frequentPatterns[i].pop(key, None)

        #Pop all non closed sequences
        for key1 in frequentPatterns[i].keys():
            for key2 in frequentPatterns[i-1].keys():
                if((key2 in key1) and (frequentPatterns[i][key1] == frequentPatterns[i-1][key2])):
                    frequentPatterns[i-1].pop(key2, None)

    allKeys = []
    allFrequentPatterns = {}
    for i in range(1, columns+1):
        for key in frequentPatterns[i].keys():
            allKeys.insert(0,key.lstrip().rstrip())
            allFrequentPatterns[key.lstrip().rstrip()] = 0
            # print key+"---"+str(frequentPatterns[i][key])

    #By this time all frequent sequences are generated
    '''*********************************************'''
    #Find the modified frequency of each sequence
    for line in blueFile:
        for key in allKeys:
            if(len(key.split(" ")) == 1):
                # print key
                key = " "+key+" "
                if(key in line):
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
                    # print line
            else:
                key = " "+key+" "
                if(key in line):
                    # key = key.lstrip().rstrip()
                    allFrequentPatterns[key.lstrip().rstrip()] = allFrequentPatterns[key.lstrip().rstrip()] + line.count(key)
                    while(key in line):
                        line = line.replace(key,' ')
    finalKeys = open("bluePatterns.txt","w")
    finalLength = 0
    for eachKey in allKeys:
        if(allFrequentPatterns[eachKey]!=0):
            finalKeys.write(eachKey.lstrip().rstrip()+"-"+str(allFrequentPatterns[eachKey])+"\n")
            finalLength = finalLength + (len(eachKey.lstrip().rstrip().split(" "))*allFrequentPatterns[eachKey])

def huffman_codes(s):
    freq = []
    i = 0
    table = []
    with open(s,"r") as s:
        for line in s:
            x = line.replace("\n",'').rstrip().split('-')
            table.append(Node([[x[0], '']], int(x[1])))
    heapq.heapify(table)
    while len(table) > 1:
        first_node = heapq.heappop(table)
        second_node = heapq.heappop(table)
        new_node = first_node.merge(second_node)
        heapq.heappush(table, new_node)
    return dict(table[0].pairs)

def huffEncode():
    s = open("redCodetable.txt","w")
    x = huffman_codes("redPatterns.txt")
    for i in x.keys():
        s.write(i+"-"+x[i]+"\n")

    s = open("greenCodetable.txt","w")
    x = huffman_codes("greenPatterns.txt")
    for i in x.keys():
        s.write(i+"-"+x[i]+"\n")

    s = open("blueCodetable.txt","w")
    x = huffman_codes("bluePatterns.txt")
    for i in x.keys():
        s.write(i+"-"+x[i]+"\n")

def huffmanEncoder(string):
    codetable = {}
    with open(string+"Codetable.txt","r") as ct:
        for line in ct:
            (pattern, code) = line.replace("\n","").split("-")
            codetable[pattern] = code

    codeTableItems = codetable.items()
    sortedCodeTable = sorted(codeTableItems,key = lambda s: len(s[0]))
    reversedSortedCodeTable = list(reversed(sortedCodeTable))
    # print reversedSortedCodeTable[0][0]
    codes = []
    for item in reversedSortedCodeTable:
        codes.append(item[0])

    # print codetable

    redEncoding = open(string+"Compressed.txt","w")
    with open(string+"ClusterEncoded.txt","r") as rc:
        for line in rc:
            line = " " + line + " "
            currentLine = line
            for key in codes:
                # print codetable[key]
                if(len(key.split(" "))==1):
                    # print "yes"
                    k1 = " "+key+" "
                    k2 = " "+key+" "
                    if((k1 in currentLine)):
                        while(k1 in currentLine):
                            currentLine = currentLine.replace(k1, " -"+codetable[key]+"- ")
                        # redEncoding.write("no-"+k1+"$$"+currentLine+"\n")
                    elif(k2 in currentLine):
                        while(k2 in currentLine):
                            currentLine = currentLine.replace(k2, " -"+codetable[key]+"- ")
                        # redEncoding.write("no-"+k2+"$$"+currentLine+"\n")
                else:
                    # print "no"
                    if(" "+key+" " in currentLine):
                        while(" "+key+" " in currentLine):
                            currentLine = currentLine.replace(" "+key+" ", " -"+codetable[key]+"- ")
                        # redEncoding.write("ues "+currentLine+"\n")
                        # print currentLprint currentLineine
            
            
            currentLine = currentLine.replace(" ","").replace("-","")
            redEncoding.write(currentLine+"\n")

def Compressor():
    for string in ['red','green','blue']:
        huffmanEncoder(string)

def huffmanDecoder(string):
    codeTable = {}
    with open(string+"Codetable.txt","r") as ct:
        for line in ct:
            (pattern, code) = line.replace("\n","").split("-")
            codeTable[code] = pattern

    redDecomp = open(string+"HuffmanDecoded.txt","w")
    with open(string+"Compressed.txt","r") as rc:
        for line in rc:
            # print line
            currentLine = ""
            cString = ""
            count = 0
            for i in line:
                # print i
                cString = cString + i
                count = 0
                for key in codeTable.keys():
                    if(cString == key):
                        count = 1
                        break
                if(count == 1):
                    # redDecomp.write(cString+"---"+codeTable[cString]+"\n")
                    # currentLine = currentLine.replace(cString,"-"+codeTable[cString]+"-")
                    currentLine= currentLine+" "+codeTable[cString]
                    # redDecomp.write(currentLine+"\n")
                    cString = ""
            redDecomp.write(currentLine.rstrip().lstrip()+"\n")
            # print currentLine
            # break

def Decoder():
    for string in ['red','green','blue']:
        huffmanDecoder(string)

def clusterDecoding(string):
    final = []
    with open(string+"HuffmanDecoded.txt","r") as s:
        m = s.read()
        for ch in m.split():
            final.append(int(ch))
    final  = np.asarray(final);
    final = final.reshape(rows,columns)
    y = blockshaped(final, blockSize, blockSize) 
    val = {}   
    for i in range(0,numberOfBlocks):
        val[i] = {}
        with open(string+"ClusterTable"+str(i)+".txt","r") as s:
            for line in s:
                mat = line.strip().split(" ")
                val[i][int(mat[0])] = int(mat[1])
    newFinal = []
    for i in range(0,numberOfBlocks):
        newY = y[i].reshape(1,blockSize*blockSize)
        yNew = []
        for element in newY[0]:
            yNew.append(val[i][int(element)])
        yNew = np.asarray(yNew)
        yNew = yNew.reshape(blockSize,blockSize)
        newFinal.append(yNew)
    newFinal = np.asarray(newFinal)
    shapedNewFinal = unblockshaped(newFinal,rows,columns)
    with open(string+"Decompressed.txt","w") as w:
        for i in range(0,rows):
            for m in shapedNewFinal[i]:
                w.write(str(m)+" ")
            w.write("\n")

def runClusterDecodingInParallel():
    proc = []
    for string in ["red","green","blue"]:
        p = Process(target=clusterDecoding, args=(string,))
        p.start()
        proc.append(p)
    for p in proc:
        p.join()

def reconstruct(rows, columns, fileName):
    data = np.zeros( (rows,columns,3), dtype=np.uint8)
    print data.shape
            
    i = 0;
    j = 0;
    with open("redDecompressed.txt","r") as rc, open("dummy.txt","w") as d:
        for line in rc:
            j = 0
            for character in line.strip().split(" "):
                data[j,i,0] = int(character)
                d.write(str(data[j,i,2])+" ")
                j = j+1;
            i = i+1
            d.write("\n")



    i = 0;
    j = 0;
    with open("greenDecompressed.txt","r") as gc, open("dummy1.txt","w") as d:
        for line in gc:
            j = 0
            for character in line.strip().split(" "):
                data[j,i,1] = int(character)
                d.write(str(data[j,i,1])+" ")
                j = j+1;
            d.write("\n")
            i = i+1



    
    i = 0;
    j = 0;
    with open("blueDecompressed.txt","r") as bc, open("dummy2.txt","w") as d:
        for line in bc:
            j = 0
            for character in line.strip().split(" "):
                data[j,i,2] = int(character)
                d.write(str(data[j,i,0])+" ")
                j = j+1;
            i = i+1
            d.write("\n")
    img = Image.fromarray(data)
    img.save(fileName+'.bmp')





res = open("results.csv","w")
res.write("BlockSize,k,Alpha,ClusteringTime,MiningTime,EncodingTime,CompressionTime,DecompressionTime,ClusterTableSize,CodeTableSize,EncodedImageSize,CompressedSize,ActualSize,JPEG Size,GIF Size,JPEG Cr,GIF Cr,Our Cr,CRP Actual,CRP JPEG,CRP GIF\n")
res.close()
for blockSize in doubling_range(32, 257):
    for numberOfClusters in range(8,33,4):
        for s in range(10,100,12):
            res = open("results.csv","a")
            image = imread("4.2.02.tiff")
            (rows,columns) = getImageSize("4.2.02.tiff")
            print str(blockSize)+"-"+str(numberOfClusters)+"-"+str(s)
            minimumSupport = (s*rows)/100.0
            numberOfBlocks = (rows*columns)/(blockSize*blockSize)
            oTS = time.clock()
            cTS = time.clock()
            runAllClusteringParallel()
            cTE = time.clock()
            print "Clustering done"
            clusteringTime = cTE-cTS
            runClusterEncodingParallel()
            print "Cluster Encoding Done"
            mTS = time.clock()
            miner()
            mTE = time.clock()
            print "Mining Done"
            mineTime = mTE-mTS
            huffEncode()
            print "Huffman Encoding Done"
            coTS = time.clock()
            Compressor()
            coTE = time.clock()
            print "Encoding Done"
            compressTime = coTE-coTS
            oTE = time.clock()
            totalCompressionTime = oTE-oTS
            print "Compression Done"
            dTS = time.clock()
            Decoder()
            runClusterDecodingInParallel()
            dTE = time.clock()
            decompressionTime = dTE-dTS
            print "Decompression Done"
            reconstruct(rows,columns,str(blockSize)+"-"+str(numberOfClusters)+"-"+str(s))
            ClTableSize = (3 * numberOfBlocks * (nob(numberOfClusters) +8) * numberOfClusters)/8000.0
            coTableSize = os.stat("redCodetable.txt").st_size+os.stat("greenCodetable.txt").st_size+os.stat("blueCodetable.txt").st_size
            coTableSize = ((coTableSize/(6*1.0)))/1000.0
            compressedImageSize = os.stat("redCompressed.txt").st_size+os.stat("greenCompressed.txt").st_size+os.stat("blueCompressed.txt").st_size
            compressedImageSize = ((compressedImageSize / 8*1.0))/1000.0 - 40

            totalSize = ClTableSize + coTableSize + compressedImageSize

            ActualSize = 786.6
            JPEGSize = 404
            GIFSize = 226
            OurCr = ActualSize/totalSize
            JPEGCr = ActualSize/JPEGSize
            GIFCr = ActualSize/GIFSize
            CRPActual = ((ActualSize-totalSize)*100)/ActualSize
            CRPJPEG = ((JPEGSize-totalSize)*100)/JPEGSize
            CRPGIF = ((GIFSize-totalSize)*100)/GIFSize
            res.write(str(blockSize)+","+str(numberOfClusters)+","+str(s)+","+str(clusteringTime)+","+str(mineTime)+","+str(compressTime)+","+str(totalCompressionTime)+","+str(decompressionTime)+","+str(ClTableSize)+","+str(coTableSize)+","+str(compressedImageSize)+","+str(totalSize)+","+str(ActualSize)+","+str(JPEGSize)+","+str(GIFSize)+","+str(JPEGCr)+","+str(GIFCr)+","+str(OurCr)+","+str(CRPActual)+","+str(CRPJPEG)+","+str(CRPGIF)+"\n")
            filelist = [ f for f in os.listdir(".") if f.endswith(".txt") ]
            for f in filelist:
                os.remove(f)
            res.close()


