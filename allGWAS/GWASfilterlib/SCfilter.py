import pandas as pd
import glob

def SCfilter(args):
    if args.screenout != None:
        screenout = glob.glob(args.screenout+r"/*.Haploview.LD.SNPblockID.screenout")
    else:
        raise Exception("Please provide the path of SNPblockID.screenout file of SNP linkage disequilibrium analysis!")

    if args.prefix != None:
        prefix = args.prefix
        Out = args.screenout + "/" + args.prefix + ".Haploview.LD.SNPblockID.screenout.SC.filterout"

    if args.anno != None:
        anno = args.anno
    else:
        raise Exception("Please provide the annotation file of all genes of the species to be analyzed!")

    if args.threshold != None:
        threshold = args.threshold


    dat = pd.read_csv(anno,"\t")
    dat = dat.fillna('')
    geneid = dat.iloc[:,10].tolist()
    sub1 = dat.iloc[:,16].tolist()
    subanno = {}
    for i in range(len(geneid)):
        if(sub1[i] == ''):
            subanno[geneid[i]] = 'Unknown'
        else:
            subanno[geneid[i]] = sub1[i].split(" //")[0].rstrip(' ')
    
    outtmp = []
    for i in screenout:
        screndat_i = pd.read_csv(i,"\t")
        screndat_i = screndat_i[screndat_i.iloc[:,4] <= threshold]
        indxi = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(screndat_i.iloc[:,0]), list(screndat_i.iloc[:,1])))
        screndat_i.index = indxi
        outtmp.append(screndat_i)
    screndat = pd.concat(outtmp)

    geneID = list(set(screndat.iloc[:,7].values.tolist()))
    filt = []
    for i in geneID:
        screndatmp = screndat[screndat.iloc[:,7] == i]
        if len(screndatmp) == 0 or len(screndatmp) < 0:
            raise Exception("There is something wrong!!!!")
        elif len(screndatmp) == 1:
            continue
        else:
            Pmin = screndatmp.iloc[:,4].min()  
            indxtmp = screndatmp.index.tolist()
            for j in range(len(indxtmp)):
                if screndatmp.iloc[j,4] > Pmin:
                    filt.append(indxtmp[j])

    screndat = screndat.drop(filt)
    Subsystem1_o = []
    for i in range(len(screndat)):
        GeneID_i = screndat.iloc[i,7]
        Subsystem1_o.append(subanno[GeneID_i])
    screndat.insert(screndat.columns.size,'Subsystem1',Subsystem1_o)
    screndat.to_csv(Out,sep="\t",header=True,index=False)



