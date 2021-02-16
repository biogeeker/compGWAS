import sys
import pickle
import logging



class writeDic:


    def __init__(self,ref_file,gpff_file,out_prefix):
        self.ref_file = ref_file
        self.gpff_file = gpff_file
        self.out_prefix = out_prefix

    

    def creatContgRefDic(self):

        big_sz=1000000000
        gene_ref_file = open(self.ref_file,'r')
        contg_ref = {}
        gene_ref0,cunt='',0

        while True:
            gene_ref1=gene_ref_file.read(big_sz)
            cunt+=1
            if gene_ref1=='':
                gene_ref=gene_ref0+'\n>'
            else:
                gene_ref=gene_ref0+gene_ref1
            st=gene_ref.find('>',0)
            while True:
                ed=gene_ref.find('\n>',st+1)
                if ed==-1:
                    gene_ref0=gene_ref[st:]
                    break
                the_blck=gene_ref[st:ed].split('\n')
                tit=the_blck[0].split()[0][1:]
                seq=''.join(the_blck[1:])
                contg_ref[tit]=seq
		
                st=ed+1
            if gene_ref1=='':
                break

        gene_ref_file.close()
        return contg_ref


    def creatAllDic(self):


        chr0,chrSt,chrEd,chrOri=1,2,3,4
        contg,contgSt,contgEd,contgOri=5,6,7,8
        gene,geneID,region,acc=9,10,11,13
        nearRNAacc,RNAacc='',{}
        contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic={},{},{},{},{},{}

        anno_file=open(self.gpff_file,'r')

        while True:
            lin=anno_file.readline().rstrip()
            if lin=='':
                break
            if lin[0] == '#':
                continue
            lst=lin.split('\t')
            theRegion=lst[region]
            theContg=lst[contg].split('.')[0]
            theGene,theGeneID=lst[gene].split('.')[0],lst[geneID]
            theAcc=lst[acc].split('.')[0]

            ST,ED,DIR=int(lst[contgSt]),int(lst[contgEd]),lst[contgOri]

            if theRegion=='CDS' and theAcc=='-':
                continue
		
            if theRegion=='GENE' and theContg not in contg_dic:
                contg_dic[theContg] =[(ST,ED,DIR,theGeneID,theGene)]  # each contig contains all genes on it.
            elif theRegion=='GENE' and theContg in contg_dic:
                contg_dic[theContg]+=[(ST,ED,DIR,theGeneID,theGene)]

            elif theRegion=='PSEUDO' and theGeneID not in mol_dic:
                mol_dic[theGeneID] =['PSEUDO']
                pseudo_dic[(theContg,theGeneID,theAcc)]=[(ST,ED,DIR,theGene)]  # This dictionary contains all pseudogenes and ranges

            elif theRegion=='PSEUDO' and theGeneID in mol_dic:
                pseudo_dic[(theContg,theGeneID,theAcc)]+=[(ST,ED,DIR,theGene)]

            elif theRegion.find('RNA')!=-1 and theGeneID not in mol_dic: 
                mol_dic[theGeneID]=['RNA']
                RNA_dic[(theContg,theGeneID,theAcc)]=[(ST,ED,DIR,theGene)]  # This dictionary contains all non-coding RNAs and ranges
            elif theRegion.find('RNA')!=-1 and theGeneID in mol_dic:
                RNA_dic[(theContg,theGeneID,theAcc)]=+[(ST,ED,DIR,theGene)]


            elif (theRegion=='UTR' or theRegion=='CDS') and theGeneID not in gene_dic: # This dictionary contains all coding genes and ranges
                gene_dic[theGeneID] =[(ST,ED,DIR,theRegion,theAcc)]
            elif (theRegion=='UTR' or theRegion=='CDS') and (ST,ED,DIR,theRegion,theAcc) not in gene_dic[theGeneID]:
                gene_dic[theGeneID]+=[(ST,ED,DIR,theRegion,theAcc)]

            if theRegion=='CDS' and (theContg,theGeneID,theAcc) not in CDS_dic:  # This dictionary contains all coding genes and coding regions
                CDS_dic[(theContg,theGeneID,theAcc)] =[(ST,ED,DIR,theGene)]
            elif theRegion=='CDS' and (theContg,theGeneID,theAcc) in CDS_dic:
                CDS_dic[(theContg,theGeneID,theAcc)]+=[(ST,ED,DIR,theGene)]

        anno_file.close()

        return contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic



    def getCDSseq(self,contg_ref,CDS_dic):

        CDSseq_dic={}       
        for CDS_kys in sorted(list(CDS_dic.keys())):
            geneSeq=''
            seqContg=CDS_kys[0]

            if seqContg not in contg_ref:
                logging.warning('the contig %s does not exist in the reference genome, ignored' % seqContg)
                break

            CDS_pos=sorted(CDS_dic[CDS_kys])
            for echPos in CDS_pos:
                geneSeq+=contg_ref[seqContg][echPos[0]-1:echPos[1]]
            CDSseq_dic[CDS_kys]=geneSeq

            if len(geneSeq)%3!=0:
                logging.warning('the protein is incomplete, recheck your annotation: ', CDS_kys)

        return CDSseq_dic



    def getPseudoSeq(self,contg_ref,pseudo_dic):

        pseudoSeq_dic={}
        for pseudo_kys in sorted(list(pseudo_dic.keys())):
            pseudoSeq=''
            pseudoContg=pseudo_kys[0]

            if pseudoContg not in contg_ref:
                logging.warning('the contig %s does not exist in the reference genome, ignored' % pseudoContg)
                break
            pseudo_pos=sorted(pseudo_dic[pseudo_kys])
            for ech_pseudoPos in pseudo_pos:
                pseudoSeq+=contg_ref[pseudoContg][ech_pseudoPos[0]-1:ech_pseudoPos[1]]
            pseudoSeq_dic[pseudo_kys]=pseudoSeq

            if len(pseudoSeq)%3==0:
                logging.warning('the pseudogene %s seems a normal coding gene,recheck your annotation. ', pseudo_kys[1] )

        return pseudoSeq_dic



    def getRNASeq(self,contg_ref,RNA_dic):

        rnaSeq_dic={}
        for RNA_kys in sorted(list(RNA_dic.keys())):
            rnaSeq=''
            rnaContg=RNA_kys[0]

            if rnaContg not in contg_ref:
                logging.warning('the contig %s does not exist in the reference genome, ignored' % rnaContg)
                break
            rna_pos=sorted(RNA_dic[RNA_kys])
            for ech_rnaPos in rna_pos:
                rnaSeq+=contg_ref[rnaContg][ech_rnaPos[0]-1:ech_rnaPos[1]]
            rnaSeq_dic[RNA_kys]=rnaSeq
        
        return rnaSeq_dic



    def dumpDic(self,contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic,CDSseq_dic,pseudoSeq_dic,rnaSeq_dic):

        out_contgDic=open(self.out_prefix +'.contg.dic','w')
        out_molDic=open(self.out_prefix +'.mol.dic','w')
        out_geneDic=open(self.out_prefix +'.gene.dic','w')
        out_CDSDic=open(self.out_prefix +'.CDS.dic','w')
        out_pseudoDic=open(self.out_prefix +'.pseudo.dic','w')
        out_rnaDic=open(self.out_prefix +'.RNA.dic','w')

        out_CDSseqDic=open(self.out_prefix +'.CDSseq.dic','w')
        out_pseudoSeqDic=open(self.out_prefix +'.pseudoSeq.dic','w')
        out_rnaSeqDic=open(self.out_prefix +'.rnaSeq.dic','w')


        #sort is not needed here!!!    
        contg_dicStr = pickle.dumps(contg_dic,0).decode()
        out_contgDic.write(contg_dicStr)

        mol_dicStr =  pickle.dumps(mol_dic,0).decode()
        out_molDic.write(mol_dicStr)

        gene_dicStr = pickle.dumps(gene_dic,0).decode()
        out_geneDic.write(gene_dicStr)

        CDS_dicStr = pickle.dumps(CDS_dic,0).decode()
        out_CDSDic.write(CDS_dicStr)

        pseudo_dicStr = pickle.dumps(pseudo_dic,0).decode()
        out_pseudoDic.write(pseudo_dicStr)

        RNA_dicStr = pickle.dumps(RNA_dic,0).decode()
        out_rnaDic.write(RNA_dicStr)

        CDSseq_dicStr = pickle.dumps(CDSseq_dic,0).decode()
        out_CDSseqDic.write(CDSseq_dicStr)

        pseudoSeq_dicStr = pickle.dumps(pseudoSeq_dic,0).decode()
        out_pseudoSeqDic.write(pseudoSeq_dicStr)

        rnaSeq_dicStr = pickle.dumps(rnaSeq_dic,0).decode()
        out_rnaSeqDic.write(rnaSeq_dicStr)


        out_contgDic.close()
        out_molDic.close()
        out_geneDic.close()
        out_CDSDic.close()
        out_pseudoDic.close()
        out_rnaDic.close()

        out_CDSseqDic.close()
        out_pseudoSeqDic.close()
        out_rnaSeqDic.close()




def parseGpff(args):
    
    if args.ref != None and args.gpff !=None and args.prefix != None:
        ref_file = args.ref
        gpff_file = args.gpff
        out_prefix = args.prefix
        logging.basicConfig(
                format='%(asctime)s \tFile \"%(filename)s" %(levelname)s: \n %(message)s \n',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=out_prefix+'.log', filemode="w", level=logging.INFO)


        allDic = writeDic(ref_file,gpff_file,out_prefix)

        contgRef = allDic.creatContgRefDic()
        contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic=allDic.creatAllDic()
        CDSseq_dic = allDic.getCDSseq(contgRef,CDS_dic)
        pseudoSeq_dic = allDic.getPseudoSeq(contgRef,pseudo_dic)
        rnaSeq_dic = allDic.getRNASeq(contgRef,RNA_dic)
        allDic.dumpDic(contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic,CDSseq_dic,pseudoSeq_dic,rnaSeq_dic)

    else:
        raise Exception("Please provide the reference genome file, gpff annotaiton file and output prefix.")




