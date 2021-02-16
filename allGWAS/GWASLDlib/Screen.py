import pandas as pd
import numpy as np
import os
import glob
import sys
import logging

def Screen(args):
    if args.allSNP != None:
        SNPfiles = glob.glob(args.allSNP+r"/*")
    else:
        raise Exception("Please provide SNP files path!")

    if args.Chi2 != None:
        Chi2 = args.Chi2
    else:
        raise Exception("Please provide the chi2 distribution file of the 2 phenotypes!")

    if args.LogisAnno != None:
        LogisAnno = args.LogisAnno
    else:
        raise Exception("Please provide the annotation file of SNP logistic regression results!")

    if args.info != None:
        INFO = args.info
    else:
        raise Exception("Please provide the Haploview.info file!")
    

    data_alt,data_ref = {},{}
    for i in range(0,len(SNPfiles)):
        i_file = SNPfiles[i]
        file_nm = i_file.split("/")[-1].split("GCA-")[-1].split(".")[0]
        file_i = pd.read_csv(i_file, "\t")
        INDEX = list(map(lambda x,y:str(x)+"_!T!E!S!T!_"+str(y), list(file_i.iloc[:,0]), list(file_i.iloc[:,1])))
        data_alt[file_nm] = pd.Series(list(file_i.iloc[:,6]), index=INDEX)
        data_alt[file_nm] = data_alt[file_nm].astype('category')
        data_alt[file_nm] = data_alt[file_nm].cat.add_categories([''])
        data_ref[file_nm] = pd.Series(list(file_i.iloc[:,4]), index=INDEX)
        data_ref[file_nm] = data_ref[file_nm].astype('category')
        data_ref[file_nm] = data_ref[file_nm].cat.add_categories([''])
    ALTData = pd.DataFrame(data_alt)
    REFData = pd.DataFrame(data_ref)
    INDEX = ALTData.index.tolist()
    ctig = list(map(lambda x:x.split("_!T!E!S!T!_")[0], INDEX)) 
    LOC = list(map(lambda x:int(x.split("_!T!E!S!T!_")[1]), INDEX))
    ALTData.insert(0, "Loc", LOC)
    ALTData.insert(0, "RefName", ctig)
    REFData.insert(0, "Loc", LOC)
    REFData.insert(0, "RefName", ctig)

    CHI2d = pd.read_csv(Chi2,sep="\t")
    for jj in INFO:
        info = jj
        filename = os.path.split(info)[1]
        contig = filename.split(args.prefix + ".")[-1].split(".Haploview.info")[0]
        LD = "".join(list(info)[0:-4]) + "LD"
        blockID = LD + ".SNPblockID"
        out = blockID + ".screenout"

        chi2flt = []
        chi2d = CHI2d[CHI2d.iloc[:,0]==contig]
        Locs = chi2d.iloc[:,1]
        freq1 = chi2d.iloc[:,2]
        freq2 = chi2d.iloc[:,3]
        freq3 = chi2d.iloc[:,4]
        freq4 = chi2d.iloc[:,5]
        for i in range(len(Locs)):
            freqs = sorted([freq1[i],freq2[i],freq3[i],freq4[i]])
            if freqs[1] < 3:
                chi2flt.append(Locs[i])

        SNPannoDic = {}
        annodat = pd.read_csv(LogisAnno,"\t")
        annodat = annodat[annodat.iloc[:,0]==contig]
        Site = annodat.iloc[:,1].values.tolist()
        for i in range(len(Site)):
            SNPannoDic[Site[i]] = annodat.iloc[i,2:9].values.tolist()

        blockdic = {}
        blockdat = pd.read_csv(blockID,sep="\t")
        blockdat = blockdat.fillna("")
        blockIDs = list(set(blockdat.iloc[:,1].values.tolist()))
        for i in blockIDs:
            blockdattmp = blockdat[blockdat.iloc[:,1]==i]
            blockdattmpNm = blockdattmp.iloc[:,0].values.tolist()
            if i == "":
                for m in blockdattmpNm:
                    if m in chi2flt:
                        continue
                    elif SNPannoDic[m][2] == "CDS_nonSynon":
                        blockdic[m] = ""
            else:
                datatmp,coeffstmp,Pvalstmp,muttmp,filtm = {},[],[],[],[]
                for mm in range(len(blockdattmpNm)):
                    m = blockdattmpNm[mm]
                    if m in chi2flt:
                        filtm.append(mm)
                    coeffstmp.append( SNPannoDic[m][0] )
                    Pvalstmp.append( SNPannoDic[m][1] )
                    muttmp.append( SNPannoDic[m][2] )
                datatmp["Loc"] = pd.Series(blockdattmpNm, index=range(len(blockdattmpNm)))
                datatmp["coeffs"] = pd.Series(coeffstmp, index=range(len(coeffstmp)))
                datatmp["muttype"] = pd.Series(muttmp, index=range(len(muttmp)))
                datatmp["Pvalues"] = pd.Series(Pvalstmp, index=range(len(Pvalstmp)))
                dataFtmp = pd.DataFrame(datatmp)
                dataFtmp = dataFtmp.drop(filtm)
                dataFtmp = dataFtmp[dataFtmp.iloc[:,2] == "CDS_nonSynon"]
                if len(dataFtmp) == 0:
                    continue
                elif len(dataFtmp) == 1:
                    key = dataFtmp.iloc[0,0]
                    blockdic[ dataFtmp.iloc[0,0] ] = i
                elif len(dataFtmp) > 1:
                    dataFtmp2 = dataFtmp[dataFtmp.iloc[:,3] == dataFtmp.iloc[:,3].min()]
                    blockdic[ dataFtmp2.iloc[0,0] ] = i
                    key = dataFtmp2.iloc[0,0]

        if len(blockdic) == 0:
            logging.basicConfig(
                    format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=blockID+'.log', filemode="w", level=logging.INFO)
            logging.warning("There is no snp screened out!")
            sys.exit()

        AltData = ALTData[ALTData.iloc[:,0]==contig]
        AltData = AltData.sort_values(by=["Loc"],ascending=True)
        AltData.index = AltData.iloc[:,1].tolist()
        AltData = AltData.iloc[:,2:]
        RefData = REFData[REFData.iloc[:,0]==contig]
        RefData = RefData.sort_values(by=["Loc"],ascending=True)
        RefData.index = RefData.iloc[:,1].tolist()
        RefData = RefData.iloc[:,2:]
        
        filtloc = list(set(AltData.index.tolist()).difference(set(blockdic.keys())))
        AltData = AltData.drop(filtloc)
        RefData = RefData.drop(filtloc)
        if args.ref != None:
            AltData[str(args.ref[0])]= np.nan
        RefData = RefData.fillna("")
        for m in range(0,len(AltData)):
            for n in range(0,RefData.shape[1]):
                if RefData.iloc[m,n] != "":
                    ref_m = RefData.iloc[m,n]
                    break
            AltData.iloc[m,:] = AltData.iloc[m,:].fillna(ref_m)

        allCunt = AltData.apply(pd.value_counts, axis=1)
        MinAle = pd.DataFrame(allCunt.idxmin(axis=1))
        MinAle.columns = ["Allele"]
        MinAle.insert(0,'Site',MinAle.index)
        MinAle.insert(0,'RefName',[contig]*len(MinAle))
        coefout,Pvalout,mutout,genout,genIDout,AAchagout,AAsitout = [],[],[],[],[],[],[]
        for i in MinAle.index.tolist():
            coefout.append( SNPannoDic[i][0] )
            Pvalout.append( SNPannoDic[i][1] )
            mutout.append( SNPannoDic[i][2] )
            genout.append( SNPannoDic[i][3] )
            genIDout.append( SNPannoDic[i][4] )
            AAchagout.append( SNPannoDic[i][5] )
            AAsitout.append( SNPannoDic[i][6] )
        MinAle.insert(MinAle.shape[1],'Allele2_Coefficient',coefout)
        MinAle.insert(MinAle.shape[1],'Allele2_Pvalue',Pvalout)
        MinAle.insert(MinAle.shape[1],'Muttype',mutout)
        MinAle.insert(MinAle.shape[1],'Gene',genout)
        MinAle.insert(MinAle.shape[1],'GeneID',genIDout)
        MinAle.insert(MinAle.shape[1],'AAchange',AAchagout)
        MinAle.insert(MinAle.shape[1],'AAsite',AAsitout)
        MinAle['AAsite'] = MinAle['AAsite'].astype(int)
        MinAle.to_csv(out,sep="\t",header=True,index=False)

