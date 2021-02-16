import pickle
import pandas as pd

def Cfilter(args):
    if args.logisre != None:
        logisre = args.logisre
        Out = logisre + ".C.filterout"
    else:
        raise Exception("Please provide the results file of CDS logistic regression!")

    if args.Chi2 != None:
        Chi2 = args.Chi2
    else:
        raise Exception("Please provide the chi2 distribution file of the 2 phenotypes!")

    if args.anno != None:
        anno = args.anno
    else:
        raise Exception("Please provide the annotation file of all genes of the species to be analyzed!")

    if args.threshold != None:
        threshold = args.threshold


    chi2flt = []
    chi2d = pd.read_csv(Chi2,sep="\t")
    Locs = chi2d.iloc[:,0]
    freq1 = chi2d.iloc[:,1]
    freq2 = chi2d.iloc[:,2]
    freq3 = chi2d.iloc[:,3]
    freq4 = chi2d.iloc[:,4]
    for i in range(len(Locs)):
        freqs = sorted([freq1[i],freq2[i],freq3[i],freq4[i]])
        if freqs[1] < 3:
            chi2flt.append(Locs[i])

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

    filtindx = []
    data = pd.read_csv(logisre,"\t")
    data = data[data.iloc[:,6] <= threshold]
    data.index = list(range(len(data)))
    for i in range(len(data)):
        if data.iloc[i,0] in chi2flt:
            filtindx.append(i)
            continue
    data = data.drop(filtindx)

    GeneID_o,Estimate_o,Pvalue_o = data.iloc[:,0].tolist(),data.iloc[:,4].tolist(),data.iloc[:,6].tolist()
    Subsystem1_o = []
    for i in range(len(GeneID_o)):
        Subsystem1_o.append(subanno[GeneID_o[i]])

    OUT = {}
    OUT["GeneID"] = pd.Series(GeneID_o, index=list(range(len(GeneID_o))))
    OUT["Allele2_Coefficient"] = pd.Series(Estimate_o, index=list(range(len(Estimate_o))))
    OUT["Allele2_Pvalue"] = pd.Series(Pvalue_o, index=list(range(len(Pvalue_o))))
    OUT["Subsystem1"] = pd.Series(Subsystem1_o, index=list(range(len(Subsystem1_o))))
    OUTF = pd.DataFrame(OUT)

    OUTF.to_csv(Out, sep="\t", header=True, index=False)






