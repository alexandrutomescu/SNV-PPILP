# Date: 21-10-'13
# Name: Karen van Rens
# Update: 13-01-14

import argparse
import subprocess
import math
from time import time, gmtime, strftime
import time
import datetime
from collections import defaultdict
import os

####################### main ##################################

qualValues = dict()

def main():
    positiondic = {}
    allMutations = []
    allQS = []
    snpQual = {}
    count = 0

    # set input parameters
    parser = argparse.ArgumentParser(description="SNV-PPILP: Refined SNV calling for tumor data using perfect phylogenies and ILP. \n Ver. 1.2")

    group = parser.add_argument_group("required arguments")
    group.add_argument("-i", dest="vcf", help="GATK's multi-sample .vcf file", required=True)
    group.add_argument("-o", dest="outputFile", help="output file", required=True)

    #parser.add_argument("-p", dest="prefix", help="name of prefix", default="")
    parser.add_argument("-ilp", dest="lpSolvePath", help="path to the directory containing lp_solve (this argument is not needed if lp_solve is in the PATH variable)", default="")
    parser.add_argument("-f", dest="f", help="an SNV non detected by GATK in one sample gets weight = (the average quality score of that SNV over all samples) - (sqrt(100/F)) * (the standard deviation of that SNV over all samples) [default=50]", default="50")
    parser.add_argument("-hc", dest="hcFlag", help="", action='store_true', default=False)


    # get input parameters
    args = parser.parse_args()
    gatkFile = args.vcf
    nr = "0" #args.prefix
    output = ""
    outputFile = args.outputFile
    lpPath = args.lpSolvePath
    hcFlag = args.hcFlag

    nSamples = 0
    sampleNames = []
    # get number of samples and their names for GATK's .vcf file
    gfile = file(gatkFile, "r")
    gline = gfile.readline()
    while gline[:6] != "#CHROM":
        gline = gfile.readline()
        continue
    gfile.close()


    cells = gline.split()
    nSamples = len(cells) - 9
    for i in range(1,nSamples+1):
        sampleNames.append(cells[8 + i])

    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Input file contains " + str(nSamples) + " samples"


    snpQual = defaultdict(list)
    # for every file; find mutations and qualities
    for idx in range(0,nSamples):    
        foundMutations, foundQual = getVCF(gatkFile, hcFlag, idx+1)
        positiondic[idx] = foundMutations
        allMutations = list(set(allMutations + foundMutations))
        allQS = list(allQS + foundQual)

        cyclicCounter = 0
        for i in range(0, len(foundMutations)):
            cyclicCounter += 1
            if cyclicCounter == 10000:
                cyclicCounter = 0
                ######## log #######
                #print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Initialized quality values for " + str(i) + " SNVs out of " + str(len(foundMutations))
            snpQual[foundMutations[i]].append(foundQual[i])
        ######## log #######
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Collected all SNVs for Sample " + str(idx+1) + ": " + str(sampleNames[idx])


    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Created binary mutation matrix"

            
    # get table through Maximum weighted independant set
    w1, c1, w2, c2, matrixdic = getMWIS(output, allMutations, positiondic, nr, lpPath, nSamples)

    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Found the columns of the mutation matrix forming a Maximum-Weight Independent Set"

    # calculate k
    k = math.sqrt(1.0/(float(args.f)/100.0))

    # calculate mean quality score
    mean = sum(allQS)/len(allQS)

    # calculate standard deviation (SD)
    difs = []
    for qs in allQS:
        difs.append((qs - mean)*(qs - mean))
    SD = math.sqrt(sum(difs)/len(allQS))

    # calculate weight of 0's
    weightZeros = mean -(k*SD)

    # add other columns to the table
    changedColumns = addToTable(output, c1, c2, w1, nr, w2, matrixdic, snpQual, weightZeros, lpPath, nSamples)

    # print changedColumns["41131586"]
   
    # print new arrangement SNP's
    sn = file(outputFile, "w")
    sn.write("Sample ID,list of SNVs (CHROM:POS)\n")

    count = 0
    tempSet1 = set(changedColumns.keys())
    for f in sampleNames:
        tempSet2 = set(positiondic[count])
        sn.write(f)
        for SNP in allMutations:
            if changedColumns:
                if SNP in tempSet1:
                    if changedColumns[SNP][1][count] == 1:
                        sn.write(", (" + SNP[0] + ":" + SNP[1] + ")")
            else:
                if SNP in tempSet2:
                    if hcFlag:
                        sn.write(", (" + SNP[0] + ":" + SNP[1] + ")")
                    else:
                        sn.write(", " + SNP)
        count += 1
        sn.write("\n")
    sn.close()


    subprocess.call("rm " + "conflicts.lp", shell=True)
    subprocess.call("rm " + "conflicts_solved.txt", shell=True)
    subprocess.call("rm " + "solveWMIS.lp", shell=True)
    subprocess.call("rm " + "solvedWMIS.txt", shell=True)

    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Job done"
    

####################### main ##################################
####################### getVCF ################################  

def getVCF(gatkFile, hcFlag, idx):	
    foundMutations = []
    foundQual = []
    
    try:
        # reads the output of GATK
        vcf = file(gatkFile, 'r')
        vcfline = vcf.readline()
        
        cyclicCounter = 0
        totalCounter = 0
        while vcfline != "":
            cyclicCounter += 1
            totalCounter += 1
            if cyclicCounter == 100000:
                cyclicCounter = 0
                ######## log #######
                # print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Processed " + str(totalCounter) + " lines of the input file for Sample " + str(idx)

            if vcfline.startswith("#") == False:
                cells = vcfline.split("\t")
                condition1 = cells[8 + idx][:3] == "1/1"
                if hcFlag:
                    condition1 = (cells[8 + idx][:3] == "1/1" or cells[8 + idx][:3] == "0/1")
                if condition1:
                    foundMutations.append((cells[0],cells[1]))
                    #qual = float(cells[8 + idx].split(":")[4].split(",")[0])
                    qual = float(cells[8 + idx].split(":")[3])
                    qualValues[(idx,(cells[0],cells[1]))] = qual
                    foundQual.append(qual)


                if cells[8 + idx][:3] == "0/0":
                    qual = 0
                    if hcFlag:
                        qual = float(cells[8 + idx].split(":")[3])
                    else:
                        qual = float(cells[8 + idx].split(":")[4].split(",")[2])                    
                    qualValues[(idx,(cells[0],cells[1]))] = qual
                if cells[8 + idx][:3] == "0/1" and False:
                    #qual = max(float(cells[8 + idx].split(":")[4].split(",")[0]), float(cells[8 + idx].split(":")[4].split(",")[2]))
                    qual = float(cells[8 + idx].split(":")[3])
                    qualValues[(idx,(cells[0],cells[1]))] = qual
    
                    #qual = float(cells[8 + idx].split(":")[4].split(",")[0])
            vcfline = vcf.readline()
        vcf.close()
    
        ######## log #######
        # print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Processed " + str(totalCounter) + " lines of the input file for Sample " + str(idx)

    except IOError:
       print "File (" + f + ") not found"

    return foundMutations, foundQual

####################### getVCF ################################
####################### getMWIS ###############################

def getMWIS(output, allMutations, positiondic, nr, lpPath, nSamples):
    matrixdic = defaultdict(list)
    w1 = []
    c1 = []
    ColumnPos = []
    selectedNodes = []
    w2 = []
    c2 = []
    countConstraints = 0


    nSamples = len(positiondic)

    positiondicEff = {}
    for i in range(0,nSamples):
        positiondicEff[i] = set(positiondic[i])

    # Get matrixdic [key=column] and [value=list of snp positions]
    cyclicCounter = 0
    totalCounter = 0
    for mutation in allMutations:
        column = []
        for sample in positiondic:
            if mutation in positiondicEff[sample]:
                column.append(1)
            else:
                column.append(0)
        
        matrixdic[tuple(column)].append(mutation)
        
        cyclicCounter += 1
        totalCounter += 1
        if cyclicCounter == 1000:
            cyclicCounter = 0
            ######## log #######
            # print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Created binary column #" + str(totalCounter)

    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Created mutation groups"

    # Create lists with weight (w1) and columns (c1)
    # weight of mutation group is sum of quality values of all entries
    for x in matrixdic:
        count = float(0)
        for mut in matrixdic[x]:
            for i in range(nSamples+1):
                if (i,mut) in qualValues:
                    count += qualValues[(i,mut)]
                    #print str((i,mut)) + str(qualValues[(i,mut)])
                    #print qualValues[(i,mut)]
                else:
                    count += 0 # improve?
        #print count
        w1.append(count)   
        c1.append(list(x))

    # #weight of mutation group is number of 1s:
    # for x in matrixdic:
    #     count = 0
    #     for i in x:
    #         if i == 1:
    #             count += 1     
    #     w1.append(len(matrixdic[x])*count)   
    #     c1.append(list(x))



    # find conflicts between matrices through ILP
    lp = file("solveWMIS.lp", 'w')
    lp.write("max: ")
    for i in range(0, len(w1)):
        lp.write(str(w1[i]) + " * x" + str(i+1))
        if i != len(w1)-1:
            lp.write(" + ")
    lp.write(";\n")
    

    for i in range(0, len(c1)):  # find all the positions of 1 in the columns and store them
        count = 0
        posOne = []
    
        for b in c1[i]:
            if b == 1:
                posOne.append(count)
            count += 1
        ColumnPos.append(set(posOne))


    noConstraints = 0
    for cp1 in ColumnPos:   # Compare positions of 1; and find conflicted matrices.
        if len(cp1) > 1:
            for cp2 in ColumnPos[ColumnPos.index(cp1)+1:]:
                if len(cp2) > 1 and len(cp1 & cp2) > 0 and len(cp1 & cp2) < len(cp1) and len(cp1 & cp2) < len(cp2):
                    countConstraints += 1
                    lp.write("x" + str(ColumnPos.index(cp1) + 1) + " + x" + str(ColumnPos.index(cp2) + 1) + " <= 1;\n")
                    noConstraints += 1


    # for i in range(len(c1)):
    #     for j in range(len(c1)):
    #         for k1 in range(nSamples):
    #             for k2 in range(nSamples):
    #                 for k3 in range(nSamples):
    #                     if c1[i][k1] == 1 and c1[j][k1] == 1 and c1[i][k2] == 1 and c1[j][k2] == 0 and c1[i][k3] == 0 and c1[j][k3] == 1:
    #                         lp.write("x" + str(i+1) + " + x" + str(j+1) + " <= 1;\n")

    
    for i in range(0, len(w1)):
        lp.write("0 <= x" + str(i+1) + " <= 1;\n")
    for i in range(0, len(w1)):
        lp.write("int x" + str(i+1) + ";\n")

    lp.close()

    # solve ILP file
    subprocess.call(lpPath + "lp_solve -s " + "solveWMIS.lp > " + "solvedWMIS.txt", shell=True)

    # read solution
    lps = file("solvedWMIS.txt","r")
    lines = lps.readlines()

    for line in lines[4:]:
        if line.endswith("1\n"):
            node = line.replace("x","").replace("1\n","").replace(" ","")
            selectedNodes.append(node)

    for sn in selectedNodes:
        w2.append(w1[int(sn)-1])
        c2.append(c1[int(sn)-1])

    return w1, c1, w2, c2, matrixdic

####################### getMWIS ###############################
####################### addToTable ############################

def addToTable(output, c1, c2, w1, nr, w2, matrixdic, snpQual, v, lpPath, nSamples): 
    matrixdicReversedConflicted = {}
    snpQualConflicted = {}
    
    # get conflicted columns
    conflictedColumns = [x for x in c1 if x not in c2]

    
    for md in matrixdic:
        if list(md) in conflictedColumns:
            for snp in matrixdic[md]:
                matrixdicReversedConflicted[snp] = md
                snpQualConflicted[snp] = snpQual[snp]
    
    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Found the columns conflicting with the MWIS"

    # if there are conflicts, change the minimum number of cells, to make all columns compatible
    if conflictedColumns != []:
        c3, w3, changedColumns = solveConflicts(output, nr, w2, c2, matrixdicReversedConflicted, snpQualConflicted, v, lpPath, nSamples)
    else:
        # if there are no conflicts
        c3 = c2
        w3 = w2
        changedColumns = {}

    return changedColumns

####################### addToTable ############################
####################### solveConflicts ########################

def solveConflicts(output, nr, w2, c2, matrixdicReversedConflicted, snpQualConflicted, v, lpPath, nSamples):
    # matrixdicReversedConflicted = { conflicted mutation : column }
    # snpQualConflicted = { conflicted mutation : list of qualityscores }
    snp = []
    sumQual = []
    averageQual = []
    columns = []
    order = []
    count = 0
    changedColumns = {}

    noConstraints = []




    #print "nSamples = " + str(nSamples)
    # sort data in "order"
    for sqc in snpQualConflicted:
        #print str(sqc)
        sumOfScores = 0
        #for s in snpQualConflicted[sqc]:
        #    sumOfScores += float(s)
        for i in range(nSamples+1):
            if (i,sqc) in qualValues:
                sumOfScores += qualValues[(i,sqc)]
        #average = sumOfScores / float(len(snpQualConflicted[sqc]))
        average = sumOfScores
        #print str(average)
        snp.append(sqc)
        sumQual.append(sumOfScores)
        averageQual.append(average)
        columns.append(list(matrixdicReversedConflicted[sqc]))

    listToSort = [(i,averageQual[i]) for i in range(len(averageQual))]
    sortedList = sorted(listToSort, key=lambda tup: tup[1])
    for i in range(len(averageQual)-1,-1,-1):
        order.append(sortedList[i][0])

    del listToSort
    del sortedList

    # for i in range(0, len(averageQual)):
    #     sumMax = -1
    #     for j in range(0, len(averageQual)):
    #         if averageQual[j] > sumMax and j not in order:
    #            sumMax =  averageQual[j]
    #            index = j
    #     order.append(index)
    #print "order = " + str(order)


    # change 0's to -1's 
    for i in range(0, len(c2)):
        for j in range(0, len(c2[0])):
            if c2[i][j] == 0:
                c2[i][j] = -1

    for col in columns: 
        for i in range(0, len(col)):
            if col[i] == 0:
                col[i] = -1



    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " " + str(len(order)) + " columns will be edited"


    cyclicCounter = 0
    totalCounter = 0
    time1 = datetime.datetime.now()

    # to make each conflicted column compatible with the matrix
    for o in order:
        columnPos = []
        count += 1
        countOnes = 0

        cyclicCounter += 1
        totalCounter += 1
        if cyclicCounter == 1000:
            ######## log #######
            cyclicCounter = 0
            time2 = datetime.datetime.now()
            delta = time2 - time1
            time1 = time2 
            # print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Edited " + str(totalCounter) + " columns \t " + "Estimated remaining time: " + str(timedelta(delta.total_seconds() * (len(order) - totalCounter)))
            print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Edited " + str(totalCounter) + " columns \t " + "Estimated time remaining: " + str(datetime.timedelta(seconds=int(delta.total_seconds() * (len(order) - totalCounter) / 1000)))


        
        # create ilp file 
        # l = file("conflicts.lp","w")
        l = ""

        # write objective function
        for i in range(0, len(columns[o])):
            if columns[o][i] == 1:
                qual = str(snpQualConflicted[snp[o]][countOnes])
                l += (qual + "*x" + str(i + 1))
                countOnes += 1 
            else: 
                l += (str(v) + "*x" + str(i + 1))
            if i == len(columns[o])-1:
                l += (";\n")
            else:
                l += (" + ")
        # for i in range(0, len(columns[o])):
        #     if (i,snp[o]) in qualValues:
        #         l.write(str(qualValues[(i,snp[o])]) + "*x" + str(i + 1))
        #         #print str(qualValues[(i,snp[o])])
        #     else:
        #         l.write(str(v) + "*x" + str(i + 1))
        #     if i == len(columns[o])-1:
        #         l.write(";\n")
        #     else:
        #         l.write(" + ")

        # write constraints
        constraints = writeConstraints(c2, columns[o])
        noConstraints.append(len(constraints))
        for cs in constraints:
            l += (cs)

        # write domain
        for i in range(0, len(columns[o])):
            l += ("0 <= x" + str(i + 1) + " <= 1;\n")
            
        # write declarations
        for i in range(0, len(columns[o])):
            l += ("int x" + str(i + 1) + ";\n")


        # create ilp file 
        lfile = file("conflicts.lp","w")
        lfile.write(l)
        lfile.close()


        # solve ilp file
        # subprocess.call(lpPath + "lp_solve -s " + "conflicts.lp > " + "conflicts_solved.txt", shell=True)
        # p = subprocess.Popen(lpPath + "lp_solve > conflicts_solved.txt <<< \"" + l + "\"", shell=True)
        p = subprocess.Popen(lpPath + "lp_solve -s " + "conflicts.lp > " + "conflicts_solved.txt",shell=True)
        p.wait()

        #exec("subprocess.call(\"" + lpPath + "lp_solve -s conflicts.lp > conflicts_solved.txt\", shell=True)")

        # read and store solution file 
        ls = file("conflicts_solved.txt","r")
        lines = ls.readlines()

        for line in lines:
            if line.endswith(" 0\n"):
                columnPos.append(line.replace("x","").replace("0\n","").replace(" ",""))


        # write original column
        changedColumns[snp[o]] = [columns[o][:]]
        
        # flip numbers according to the ilp solution        
        for b in columnPos:
            b = int(b) - 1
            if columns[o][b] == -1:
                columns[o][b] = 1
            elif columns[o][b] == 1:
                columns[o][b] = -1
        
        # add column to matrix
        if columns[o] in c2:
            w2[c2.index(columns[o])] = w2[c2.index(columns[o])] + 1
        else:
            c2.append(columns[o])
            w2.append(1)

        changedColumns[snp[o]] = changedColumns[snp[o]] + [columns[o]]

    ######## log #######
    print strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " Corrected " + str(totalCounter) + " columns"

    # change -1's to 0's
    for j in range(0, len(c2)):
        for k in range(0, len(c2[0])):
            if c2[j][k] == -1:
                c2[j][k] = 0

    
    for chcol in changedColumns:
        for j in range(0, len(changedColumns[chcol])):
            for k in range(0, len(changedColumns[chcol][j])):
                if changedColumns[chcol][j][k] == -1:
                    changedColumns[chcol][j][k] = 0

    #print "avg noConstraints in columns = " + str(sum(noConstraints) / len(noConstraints))

    return c2, w2, changedColumns

####################### solveConflicts ########################
####################### writeConstraints ######################

def writeConstraints(table, cc):
    constraints = []
    # for every combination of three cells per row      
    for column in table:
        indexFirstRow = 0
        for row1 in column:
            indexFirstRow += 1
            indexSecondRow = indexFirstRow
            for row2 in column[indexFirstRow:]:
                indexSecondRow += 1
                indexThirdRow = indexSecondRow
                for row3 in column[indexSecondRow:]:
                    indexThirdRow += 1


                    # if there is a conflict
                    if row1 + row2 + row3 == 1:
                        
                        # find indexes of rows, a and b have a row with value 1, c has a row with value 0
                        if row1 == -1:
                            a = indexSecondRow
                            b = indexThirdRow
                            c = indexFirstRow
                        elif row2 == -1:
                            a = indexFirstRow
                            b = indexThirdRow
                            c = indexSecondRow
                        elif row3 == -1:
                            a = indexFirstRow
                            b = indexSecondRow
                            c = indexThirdRow

                        # write contstraints
                        constraint = str(2*cc[a - 1]) + " * x" + str(a) + " - " + str(cc[a - 1]) + " + " + str(-2*cc[b - 1]) + " * x" + str(b) + " + " + str(cc[b - 1]) + " + " + str(2*cc[c - 1]) + " * x" + str(c) + " - " + str(cc[c - 1]) + " <= 2;\n"
                        constraint2 = str(-2*cc[a - 1]) + " * x" + str(a) + " + " + str(cc[a - 1]) + " + " + str(2*cc[b - 1]) + " * x" + str(b) + " - " + str(cc[b - 1]) + " + " + str(2*cc[c - 1]) + " * x" + str(c) + " - " + str(cc[c - 1]) + " <= 2;\n"
                        if constraint not in constraints:
                            constraints.append(constraint)
                        if constraint2 not in constraints:
                            constraints.append(constraint2)


    return constraints
                        
####################### writeConstraints ######################

main()


