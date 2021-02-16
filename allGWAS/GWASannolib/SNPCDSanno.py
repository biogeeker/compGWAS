import pandas as pd
import numpy as np

def SNPCDSanno(args):

    if args.logisre != None:
        logisre = args.logisre
        out_file = logisre + ".annotation"
    else:
        raise Exception("Please provide the results file of SNP or CDS logistic regression!")

    if args.anno != None:
        anno = args.anno
    else:
        raise Exception("Please provide the annotation file of all SNPs or all genes of the species to be analyzed!")

    if args.type == "SNP":
        indx = 0
        indx2 = 1
        olst = [10,11,12,14,15]
        Olst = [5,6]
        nlst = ["MutType","Gene","GeneID","AAchange","AAsite"]
        Nlst = ["Allele2_Coefficient","Allele2_Pvalue"]
    elif args.type == "CDS":
        indx = 10
        olst = [9,14,15,16,17]
        Olst = [4,6]
        nlst = ["Gene","Function","Subsystem1","Subsystem2","Subsystem3"]
        Nlst = ["Allele2_Coefficient","Allele2_Pvalue"]
    else:
        raise Exception("Please provide the type of annotation, SNP or CDS!")

    
    inputf = open(anno,"r")
    data = inputf.readlines()
    annodic = {}
    for i in range(0,len(data)):
        res = []
        cols = data[i].split("\n")[0].split("\t")
        if args.type == "CDS":
            key = str(cols[indx])
        elif args.type == "SNP":
            key = str(cols[indx]) + "_!T!E!S!T!_" + str(cols[indx2])
        else:
            raise Exception("Please provide the type of annotation, SNP or CDS!")
        for m in olst:
            res.append(cols[m])
        annodic[key] = res
    inputf.close()

    data = pd.read_csv(logisre,"\t")
    if len(data) == 0:
        raise Exception("There is no SNP or CDS to annotate!")

    if Olst != None and Nlst != None:
        logisdic = {}
        for i in range(len(data)):
            res = []
            if args.type == "CDS":
                key = str(data.iloc[i,0])
            elif args.type == "SNP":
                key = str(data.iloc[i,0]) + "_!T!E!S!T!_" + str(data.iloc[i,1])
            for m in Olst:
                res.append(data.iloc[i,m])
            logisdic[key] = res
    DATA,result = {},[]
    for i in range(len(data)):
        if args.type == "CDS":
            key = data.iloc[i,0]
        elif args.type == "SNP":
            key = str(data.iloc[i,0]) + "_!T!E!S!T!_" + str(data.iloc[i,1])
        if isinstance(key, str) == False:
            key = str(key)
        if Olst != None and Nlst != None:
            DATA[key] = pd.Series(logisdic[key]+annodic[key], index=Nlst+nlst)
    dataF=pd.DataFrame(DATA)
    out = pd.DataFrame(dataF.values.T, index=dataF.columns, columns=dataF.index)
    if args.type == "CDS":
        out.insert(0,"Site",out.index.tolist())
    elif args.type == "SNP":
        INDEX = out.index.tolist()
        ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], INDEX))
        loc = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), INDEX))
        out.insert(0,"Site",loc)
        out.insert(0,"RefName",ctig)
    out.to_csv(out_file,sep="\t",header=True,index=False)



