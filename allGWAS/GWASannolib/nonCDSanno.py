import pickle
import pandas as pd
import re
import logging


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()


def calculate_overlap(window_start,window_end,nonCDS_start,nonCDS_end):
    a1 = int(window_start)
    a2 = int(window_end)
    b1 = int(nonCDS_start)
    b2 = int(nonCDS_end)
    if(b1<=a1 and a2 <= b2):
        overlap = a2 - a1 + 1
    if(b1>=a1 and a2 >= b2):
        overlap = b2 - b1 + 1
    if ( (b1 >= a1 and b1 <= a2 and b2 > a2) ):
        overlap = a2 - b1 + 1
    if( (b2 >= a1 and b2 <= a2 and b1 < a1) ):
        overlap = b2 - a1 + 1
    return overlap


def nonCDSanno(args):
    if args.CDSdic != None:
        CDSdic = args.CDSdic
    else:
        raise Exception("Please provide CDS annotation dictionaries of the species to be analyzed!")

    if args.genome != None:
        genome = args.genome
    else:
        raise Exception("Please provide the FASTA format genome file of the species to be analyzed!")

    if args.size != None:
        size = int(args.size)
    else:
        raise Exception("Please provide the size of the window!")

    if args.logisre != None:
        logisre = pd.read_csv(args.logisre,"\t")
        out_file = open(args.logisre+".annotation","w")
    else:
        raise Exception("Please provide the results file of nonCDS logistic regression!")

    logging.basicConfig(
            format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
            datefmt='%a, %d %b %Y %H:%M:%S',
            filename=args.logisre+'.annotation.log', filemode="w", level=logging.INFO)


    anno_dic = {}
    for m in range(0,len(CDSdic)):
        anno_m = open(CDSdic[m],"r")
        anno_m_dic = pickle.load(StrToBytes(anno_m))
        anno_dic.update(anno_m_dic)
        anno_m.close()
    contigs = []
    for i in anno_dic.keys():
        contigs.append(i[0])
    contigs = list(set(contigs))

    genome_file = open(genome,"r")
    genome_tmp = genome_file.read()
    CtigLendic = {}
    for i in contigs:
        tag = i
        re_genome = re.compile(">"+tag+" "+ '.+\n' + '((\w+\n)+)')
        genome_new = list(re_genome.findall(genome_tmp)[0][0])
        while "\n" in genome_new:
            genome_new.remove("\n")
        CtigLendic[tag] = len(genome_new)

    result = 'RefName\tSite\tSub\tOverlap\tMutType\tGene\tGeneID\tAllele2_Coefficient\tAllele2_Pvalue\n'
    for jj in contigs:
        Loc,Sub,Estimate_dic,Pvalue_dic,OUT = [],[],{},{},[]
        DIP_lines = logisre[logisre.iloc[:,0]==jj]
        for i in range(0,len(DIP_lines)):
            loc = int(DIP_lines.iloc[i,1])
            sub = int(DIP_lines.iloc[i,2])
            Estimate = DIP_lines.iloc[i,7]
            Pvalue = DIP_lines.iloc[i,6]
            Loc.append(int(loc))
            Sub.append(int(sub))
            Estimate_dic[(loc,sub)] = Estimate
            Pvalue_dic[(loc,sub)] = Pvalue

        newdic = {}
        for i in anno_dic.keys():
            if i[0] == jj:
                start_i = int(anno_dic[i][0][0])
                newdic[start_i] = i
        keys = sorted(list(newdic.keys()))

        # filter the start site of gene whose both start site and end site are in an another gene
        # creat a dic (anno_dic_2) which contains all CDS regions
        anno_dic_2,filter_start = {},[]
        for mm in keys:
            i = newdic[mm]
            start_i = int(anno_dic[i][0][0])
            end_i = int(anno_dic[i][0][1])
            geneID_left_i = i[1]
            geneID_right_i = i[1]
            gene_left_i = anno_dic[i][0][3]
            gene_right_i = anno_dic[i][0][3]
            if(anno_dic[i][0][2] == "+"):
                ori_left_i = "5"
                ori_right_i = "3"
            if (anno_dic[i][0][2] == "-"):
                ori_left_i = "3"
                ori_right_i = "5"
            for n in keys:
                m = newdic[n]
                start_m = int(anno_dic[m][0][0])
                end_m = int(anno_dic[m][0][1])
                geneID_left_m = m[1]
                geneID_right_m = m[1]
                gene_left_m = anno_dic[m][0][3]
                gene_right_m = anno_dic[m][0][3]
                if (anno_dic[m][0][2] == "+"):
                    ori_left_m = "5"
                    ori_right_m = "3"
                if (anno_dic[m][0][2] == "-"):
                    ori_left_m = "3"
                    ori_right_m = "5"
                if (start_i < start_m and end_i > end_m):
                    filter_start.append(int(start_m))
                if (start_i < start_m and end_i < end_m and end_i >= start_m):
                    filter_start.append(int(start_m))
                    end_i = end_m
                    geneID_right_i = geneID_right_m
                    gene_right_i = gene_right_m
                    ori_right_i = ori_right_m
                if (start_i > start_m and start_i <= end_m and end_i >= end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                    geneID_left_i = geneID_left_m
                    gene_left_i = gene_left_m
                    ori_left_i = ori_left_m
                if (start_i > start_m and end_i < end_m):
                    filter_start.append(int(start_i))
                    start_i = start_m
                    geneID_left_i = geneID_left_m
                    gene_left_i = gene_left_m
                    ori_left_i = ori_left_m
                    end_i = end_m
                    geneID_right_i = geneID_right_m
                    gene_right_i = gene_right_m
                    ori_right_i = ori_right_m

                if (start_i == start_m and end_i == end_m and anno_dic[m] != anno_dic[i]):
                    logging.info('There is an another gene {1}, which possesses the same start and end site with the gene {0}.'.format(i[1],m[1]))
                if (start_i == start_m and end_i < end_m and end_i >= start_m):
                    logging.info('There is an another gene {1}, whose start site is same but end site is bigger than the gene {0}.'.format(i[1],m[1]))
                if (start_i > start_m and start_i <= end_m and end_i == end_m):
                    logging.info('There is an another gene {1}, whose end site is same but start site is smaller than the gene {0}.'.format(i[1],m[1]))
            if (start_i not in filter_start):
                anno_dic_2[start_i] = [start_i,end_i,geneID_left_i,geneID_right_i,gene_left_i,gene_right_i,ori_left_i,ori_right_i]

        dic2_keys = sorted(list(anno_dic_2.keys()))

        # special_region contains all windows which are tatally in CDS (pseudogene or RNA regions)
        special_region = {}
        for n in range(len(Loc)):
            i = Loc[n]
            sub_i = Sub[n]
            for m in dic2_keys:
                if(int(i)>=anno_dic_2[m][0] and int(int(i)+size-1) <=anno_dic_2[m][1]):
                    logging.info('The window {0} is totally in CDS region.\n {1}'.format(i,anno_dic_2[m]))
                    special_region[(i,sub_i)] = [str(i),str(sub_i),str(size),"inter-gene-"+anno_dic_2[m][6]+"-"+anno_dic_2[m][7],anno_dic_2[m][4]+"||"+anno_dic_2[m][5],anno_dic_2[m][2]+"||"+anno_dic_2[m][3],str(Estimate_dic[(i,sub_i)]),str(Pvalue_dic[(i,sub_i)])]

        # creat a dic (nonCDS_anno_dic) which contains all nonCDS regions
        nonCDS_anno_dic = {}
        nonCDS_anno_dic[1] = [1,dic2_keys[0]-1,'',anno_dic_2[dic2_keys[0]][2],'',anno_dic_2[dic2_keys[0]][4],'',anno_dic_2[dic2_keys[0]][6]]
        for i in range(0,len(dic2_keys)-1):
            left = anno_dic_2[dic2_keys[i]][1]+1
            geneID_left = anno_dic_2[dic2_keys[i]][3]
            gene_left = anno_dic_2[dic2_keys[i]][5]
            ori_left = anno_dic_2[dic2_keys[i]][7]
            right = anno_dic_2[dic2_keys[i+1]][0]-1
            geneID_right = anno_dic_2[dic2_keys[i+1]][2]
            gene_right = anno_dic_2[dic2_keys[i+1]][4]
            ori_right = anno_dic_2[dic2_keys[i+1]][6]
            nonCDS_anno_dic[left] = [left,right,geneID_left,geneID_right,gene_left,gene_right,ori_left,ori_right]
        if (anno_dic_2[dic2_keys[-1]][1] < CtigLendic[jj]):
            nonCDS_anno_dic[anno_dic_2[dic2_keys[-1]][1]+1] = [anno_dic_2[dic2_keys[-1]][1]+1,CtigLendic[jj],anno_dic_2[dic2_keys[-1]][3],'',anno_dic_2[dic2_keys[-1]][5],'',anno_dic_2[dic2_keys[-1]][7],'']
    
        del anno_dic_2
        del anno_dic

        # annotate all windows
        for n in range(len(Loc)):
            i,sub_i,count = Loc[n],Sub[n],0
            out_i,over,gene,geneID,ori = [],[],[],[],["inter-gene-"]
            if ( (i,sub_i) in special_region.keys() ):
                out_i = special_region[(i,sub_i)]
                OUT.append("\t".join(out_i))
                continue
            for m in nonCDS_anno_dic.keys():
                if(nonCDS_anno_dic[m][0]<int(i) and nonCDS_anno_dic[m][1]>int(int(i)+size-1)):
                    overlap = calculate_overlap(int(i), int(int(i)+size-1), nonCDS_anno_dic[m][0], nonCDS_anno_dic[m][1])
                    over.append(str(overlap))
                    ori.append(nonCDS_anno_dic[m][6] + "-" + nonCDS_anno_dic[m][7])
                    gene.append(nonCDS_anno_dic[m][4])
                    gene.append(nonCDS_anno_dic[m][5])
                    geneID.append(nonCDS_anno_dic[m][2])
                    geneID.append(nonCDS_anno_dic[m][3])
                    continue
                if ((count == 0) and ((nonCDS_anno_dic[m][0] >= int(i) and nonCDS_anno_dic[m][0] <= int(int(i) + size-1)) or (nonCDS_anno_dic[m][1] >= int(i) and nonCDS_anno_dic[m][1] <= int(int(i) + size-1)) )):
                    count =count+1
                    overlap = calculate_overlap(int(i), int(int(i)+size-1), nonCDS_anno_dic[m][0], nonCDS_anno_dic[m][1])
                    over.append(str(overlap))
                    ori.append(nonCDS_anno_dic[m][6]+"-"+nonCDS_anno_dic[m][7])
                    gene.append(nonCDS_anno_dic[m][4])
                    gene.append(nonCDS_anno_dic[m][5])
                    geneID.append(nonCDS_anno_dic[m][2])
                    geneID.append(nonCDS_anno_dic[m][3])
                    continue
                if ((count > 0) and ((nonCDS_anno_dic[m][0] >= int(i) and nonCDS_anno_dic[m][0] <= int(int(i) + size-1)) or (nonCDS_anno_dic[m][1] >= int(i) and nonCDS_anno_dic[m][1] <= int(int(i) + size-1)) )):
                    count = count + 1
                    overlap = calculate_overlap(int(i), int(int(i)+size-1), nonCDS_anno_dic[m][0], nonCDS_anno_dic[m][1])
                    over.append(str(overlap))
                    ori.append("T"+nonCDS_anno_dic[m][6]+"-"+nonCDS_anno_dic[m][7])
                    gene.append(nonCDS_anno_dic[m][5])
                    geneID.append(nonCDS_anno_dic[m][3])
                    continue
            out_i.append(jj)
            out_i.append(str(i))
            out_i.append(str(sub_i))
            out_i.append("|".join(over))
            out_i.append("".join(ori))
            out_i.append("||".join(gene))
            out_i.append("||".join(geneID))
            out_i.append(str(Estimate_dic[(i,sub_i)]))
            out_i.append(str(Pvalue_dic[(i,sub_i)]))
            OUT.append("\t".join(out_i))
        result = result + "\n".join(OUT)
    out_file.write(result)
    out_file.close()




