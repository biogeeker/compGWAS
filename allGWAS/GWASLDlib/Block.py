import pandas as pd


def Block(args):
    if args.info != None:
        INFO = args.info
    else:
        raise Exception("Please provide the Haploview.info file!")

    if args.threshold2 != None:
        threshold = args.threshold2

    if args.maxdist != None:
        maxdist = args.maxdist

    for jj in INFO:
        info = jj
        LD = "".join(list(info)[0:-4]) + "LD"
        out1 = LD + ".blocks"
        out2 = LD + ".SNPblockID"

        infodat = pd.read_csv(info,sep="\t",header=None)
        infodat.columns = ["name","loc"]
        infodat.sort_values("loc",inplace=True)
        infodat["name"] = infodat["name"].astype(str)
        infodic = {}
        index = infodat.iloc[:,0].values.tolist()
        loc = infodat.iloc[:,1].values.tolist()
        for i in range(len(index)):
            infodic[index[i]]=loc[i]

        LDdic = {}
        LDdat = pd.read_csv(LD,"\t")
        L1 = LDdat.iloc[:,0].values.tolist()
        L2 = LDdat.iloc[:,1].values.tolist()
        cutoffCol = LDdat.iloc[:,4].values.tolist()
        for i in range(len(L1)):
            LDdic[(str(L1[i]),str(L2[i]))] = cutoffCol[i]
            LDdic[(str(L2[i]),str(L1[i]))] = cutoffCol[i]
        del LDdat

        blocks,locs = [],[]
        infonames = infodat.iloc[:,0].values.tolist()
        tag,num,account = 0,0,0
        for i in range(0,len(infonames)-1):
            key = (infonames[i],infonames[i+1])
            if key in LDdic.keys():
                if( (LDdic[key] >= float(threshold)) and ((infodic[infonames[i+1]] - infodic[infonames[tag]] + 1) <= (float(maxdist)*1000)) ):
                    account = account + 1
                else:
                    if account > 0:
                        num = num + 1
                        blocks.append("block"+str(num))
                        locs.append(",".join(infonames[tag:i+1]))
                        tag = i + 1
                    else:
                        tag = tag + 1
                    account = 0 
            else:
                if account > 0:
                    num = num + 1
                    blocks.append("block"+str(num))
                    locs.append(",".join(infonames[tag:i+1]))
                    tag = i + 1
                else:
                    tag = tag + 1
                account = 0

        data = {}
        data["BlockID"] = pd.Series(blocks,index=range(len(blocks)))
        data["Names"] = pd.Series(locs,index=range(len(locs)))
        dataF = pd.DataFrame(data)
        dataF.to_csv(out1,sep="\t",header=True,index=False)

        locblockdic = {}
        for i in range(len(locs)):
            locstmp = locs[i].split(",")
            for m in locstmp:
                locblockdic[m] = blocks[i]

        snploc,blockID = [],[]
        for i in infonames:
            snploc.append( infodic[i] )
            if i in locblockdic.keys():
                blockID.append( locblockdic[i] )
            else:
                blockID.append("")
        data = {}
        data["Site"] = pd.Series(snploc,index=range(len(snploc)))
        data["BlockID"] = pd.Series(blockID,index=range(len(blockID)))
        dataF = pd.DataFrame(data)
        dataF.to_csv(out2,sep="\t",header=True,index=False)

