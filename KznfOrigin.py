import pandas as pd
import numpy as np
from Bio import pairwise2
from selenium import webdriver

'''
Using the similarity of KZN Fs domain sequence, the occurrence time sequence of KZFPs, 
and the location information of KZFPs on the human genome, the KRAB domain and C2H2 domain of KZFPs are traced
'''

def get_kznfDomainInfor():
    '''
    Use the comment information of Uniprot to obtain the location of each KZFPs structure domain on Sequence and store it
    :return: id_domainDone_dic.npy
    '''
    KznfGeneAge_dic = np.load(r'C:\Users\wu\Desktop\Evolution_paper\KznfGeneAge_dic.npy',allow_pickle=True)[()]
    EvolutionAge_dic = KznfGeneAge_dic['Evolution']
    KeyAA_dic = KznfGeneAge_dic['KeyAA']
    geneId_dic = KznfGeneAge_dic['geneId']
    remove_list = []
    for kznf in EvolutionAge_dic.keys():
        if kznf not in geneId_dic.keys():
            remove_list.append(kznf)
    for i in remove_list:
        EvolutionAge_dic.pop(i)
    EvoKznf_list = list(geneId_dic.keys())
    option = webdriver.ChromeOptions()
    prefs = {"profile.default_content_setting_values.automatic_downloads": 1,
             'profile.default_content_settings.popups': 0,
             'download.default_directory': r'C:\Users\wu\Desktop\506_new\analysis\ensemble_orthology',
             "safebrowsing.enabled": True}
    option.add_experimental_option('prefs', prefs)
    # option.add_argument('headless')
    browser = webdriver.Chrome(executable_path=r'C:\Users\wu\AppData\Local\Google\Chrome\Application\chromedriver.exe',
                               options=option)
    browser.maximize_window()
    id_domain_dic = {}
    for kznf in EvoKznf_list:
        Uniprot_id = geneId_dic[kznf]
        try:
            full_url = f'https://www.uniprot.org/uniprot/{Uniprot_id}.txt'
            # full_url = f'https://www.uniprot.org/uniprot/G3MXI2.txt'
            browser.get(full_url)
            txt_file = browser.find_element_by_xpath('/html/body/pre').text
            txt_list = txt_file.split('\n')
            infro_dic = {}
            for idx, line in enumerate(txt_list):
                if line[:2] == 'ID':
                    seqLen = line.split(';         ')[-1].strip()
                    infro_dic['seqLen'] = seqLen
                if 'FT                   /note="KRAB"' in line or 'FT                   /note="KRAB 1"' in line \
                        or 'FT                   /note="KRAB-related' in line:
                    krab_record = txt_list[idx - 1][21:].strip()
                    infro_dic[krab_record] = 'KRAB'
                if '/note="C2H2-type' in line and 'degenerate' not in line:
                    c2h2_record = txt_list[idx - 1][21:].strip()
                    infro_dic[c2h2_record] = 'C2h2'
                if '/note="C2H2-type' in line and 'degenerate' in line:
                    deC2h2_record = txt_list[idx - 1][21:].strip()
                    infro_dic[deC2h2_record] = 'dec2h2'
                if '/note="SCAN' in line:
                    Scan_record = txt_list[idx - 1][21:].strip()
                    infro_dic[Scan_record] = 'scan'
                if line[:2] == 'SQ':
                    seq = ''.join(txt_list[idx+1:-1])
                    seq = seq.replace(' ','')
                    infro_dic['Protein_seq'] = seq
            id_domain_dic[Uniprot_id] = infro_dic
        except:
            print(Uniprot_id)
    id_domainDone_dic = {}
    for id in id_domain_dic.keys():
        singleKznf_domain = id_domain_dic[id]
        if 'Reviewed' in singleKznf_domain['seqLen']:
            singleKznf_domain['seqLen'] = singleKznf_domain['seqLen'].split(';        ')[1]
        id_domainDone_dic[id] = singleKznf_domain
    id_domainSeq_dic = {}
    for id in id_domainDone_dic.keys():
        singe_domainSeq_dic = {}
        id_inforDomain_dic = id_domainDone_dic[id]
        kznf_seq = id_inforDomain_dic['Protein_seq']
        C2h2_list = []
        dec2h2_list = []
        for i in id_inforDomain_dic.keys():
            if id_inforDomain_dic[i] == 'KRAB':
                pos1 = int(i.split('..')[0])
                pos2 = int(i.split('..')[1])
                Krab_seq = kznf_seq[pos1-1:pos2]
                singe_domainSeq_dic['krab'] = Krab_seq
            if  id_inforDomain_dic[i] == 'C2h2':
                pos1 = int(i.split('..')[0])
                pos2 = int(i.split('..')[1])
                C2h2_seq = kznf_seq[pos1-1:pos2]
                C2h2_list.append(C2h2_seq)
                singe_domainSeq_dic[f'c2h2_{len(C2h2_list)}'] = C2h2_seq
            if  id_inforDomain_dic[i] == 'dec2h2':
                pos1 = int(i.split('..')[0])
                pos2 = int(i.split('..')[1])
                deC2h2_seq = kznf_seq[pos1-1:pos2]
                dec2h2_list.append(deC2h2_seq)
                singe_domainSeq_dic[f'dec2h2_{len(dec2h2_list)}'] = deC2h2_seq
        id_domainSeq_dic[id] = singe_domainSeq_dic
    np.save(r'C:\Users\wu\Desktop\Evolution_paper\id_domainDone_dic.npy',id_domainDone_dic)

def AlignOrigin(id,aimSeq,idSeq_dic):
    '''
    For the input sequence, take the result with consistency ranking in the top 40
    :return: sorted_SeqSimility_list[:40]
    '''
    SeqID_dic = {}
    for key,value in idSeq_dic.items():
        SeqID_dic[value] = key
    Seqlist = SeqID_dic.keys()
    IdSimility_dic = {}
    def ProteinSimilarity(seq1,
                          seq2):
            score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
            # print(global_align[0])
            percent_match = (score / len(seq1)) * 100
            return percent_match
    for i in Seqlist:
        Similarity = ProteinSimilarity(aimSeq,i)
        # print(id,SeqID_dic[i])
        if id not in SeqID_dic[i]:
            IdSimility_dic[SeqID_dic[i]] = Similarity
    sorted_SeqSimility_list = sorted(IdSimility_dic.items(), key=lambda d: d[1], reverse=True)
    for idx,item in enumerate(sorted_SeqSimility_list):
        sorted_SeqSimility_list[idx] = (item[1],item[0])
    return sorted_SeqSimility_list[:40]

def KznfDomainOrigin():
    '''
    Summarize and process the results of structure domain comparison
    :return: KznfDomainOrigin_dic
    '''
    id_domainSeq_dic = np.load(r'\data\id_domainSeq_dic.npy',allow_pickle=True)[()]
    KznfGeneAge_dic = np.load(r'\data\KznfGeneAge_dic.npy',allow_pickle=True)[()]
    EvolutionAge_dic = KznfGeneAge_dic['Evolution']
    geneId_dic = KznfGeneAge_dic['geneId']
    remove_list = []
    for kznf in EvolutionAge_dic.keys():
        if kznf not in geneId_dic.keys():
            remove_list.append(kznf)
    for i in remove_list:
        EvolutionAge_dic.pop(i)
    allKznf_list = [x for x in EvolutionAge_dic.keys()]
    id_krabSeq_dic = {}
    for id in id_domainSeq_dic.keys():
        id_krabSeq_dic[id] = id_domainSeq_dic[id]['krab']
    id_c2h2_dic = {}
    for id in id_domainSeq_dic.keys():
        singleKznf_dic = id_domainSeq_dic[id]
        c2h2_list = []
        for domain in singleKznf_dic.keys():
            if domain.split('_')[0] == 'c2h2':
                c2h2_list.append(singleKznf_dic[domain])
                id_c2h2_dic[f'{id}_c2h2_{len(c2h2_list)}'] = singleKznf_dic[domain]
    id_dec2h2_dic = {}
    for id in id_domainSeq_dic.keys():
        singleKznf_dic = id_domainSeq_dic[id]
        dec2h2_list = []
        for domain in singleKznf_dic.keys():
            if domain.split('_')[0] == 'dec2h2':
                dec2h2_list.append(singleKznf_dic[domain])
                id_dec2h2_dic[f'{id}_dec2h2_{len(dec2h2_list)}'] = singleKznf_dic[domain]
    id_allC2h2_dic = dict(id_c2h2_dic,**id_dec2h2_dic)
    # Compare similarities and determine possible sources
    KznfDomainOrigin_dic = {}
    for kznf in allKznf_list:
        try:
            singlekznf_allDomain_simi_dic = {}
            kznfDomain_dic = id_domainSeq_dic[geneId_dic[kznf]]
            for domain in kznfDomain_dic.keys():
                # print(domain)
                if domain == 'krab':
                    krab_seq = kznfDomain_dic[domain]
                    Similiarity_list = AlignOrigin(geneId_dic[kznf],krab_seq,id_krabSeq_dic)
                    singlekznf_allDomain_simi_dic[domain] = Similiarity_list
                if domain.split('_')[0] == 'c2h2':
                    c2h2_seq = kznfDomain_dic[domain]
                    Similiarity_list = AlignOrigin(geneId_dic[kznf], c2h2_seq, id_allC2h2_dic)
                    singlekznf_allDomain_simi_dic[domain] = Similiarity_list
                if domain.split('_')[0] == 'dec2h2':
                    dec2h2_seq = kznfDomain_dic[domain]
                    Similiarity_list = AlignOrigin(geneId_dic[kznf], dec2h2_seq, id_allC2h2_dic)
                    singlekznf_allDomain_simi_dic[domain] = Similiarity_list
            KznfDomainOrigin_dic[kznf] = singlekznf_allDomain_simi_dic
            print(kznf)
        except:
            print('failed',kznf)
    np.save(r'\data\KznfDomainOriginFromAllC2h2_AllKznf_dic.npy',KznfDomainOrigin_dic)

def geneDistence():
    '''
    Processing gene distance information between KZFPs
    :return: kznfDistance.csv
    '''
    id_domainSeq_dic = np.load(r'\data\id_domainDone_dic.npy',allow_pickle=True)[()]
    KznfGeneAge_dic = np.load(r'\data\KznfGeneAge_dic.npy',allow_pickle=True)[()]
    EvolutionAge_dic = KznfGeneAge_dic['Evolution']
    geneId_dic = KznfGeneAge_dic['geneId']
    Idgene_dic = {}
    for keys,value in geneId_dic.items():
        Idgene_dic[value] = keys
    kznfGene_list = list(geneId_dic.keys())
    geneDisAnno_df = pd.read_csv(r'\data\hg38.refGene.gtf',sep='\t').iloc[:,[0,3,-1]]
    geneDisAnno_df['gene'] = geneDisAnno_df.iloc[:,-1].apply(lambda x:x.split('gene_name "')[1][:-2])
    geneDisAnno_df.drop_duplicates(subset=['gene'], keep='first', inplace=True)
    kznfDisAnno_df = geneDisAnno_df.loc[geneDisAnno_df['gene'].isin(kznfGene_list)]
    kznfDisAnno_df.to_csv(r'\data\kznfDistance.csv')

def getKznfDistanceScore(kznfDisAnno_df,kznf1,kznf2):
    '''
    Convert the distance between KZFPs into a distance fraction
    :return: KznfDistanceScore
    '''
    kznf1_df = kznfDisAnno_df.loc[kznfDisAnno_df['gene'] == kznf1]
    kznf1_chr,kznf1_pos= str(kznf1_df['chr'].values[0]),int(kznf1_df['pos'].values)
    kznf2_df = kznfDisAnno_df.loc[kznfDisAnno_df['gene'] == kznf2]
    kznf2_chr, kznf2_pos = str(kznf2_df['chr'].values[0]), int(kznf2_df['pos'].values)
    if kznf1_chr == kznf2_chr:
        if abs(kznf1_pos - kznf2_pos) < 4000000:
            score  = 4 - abs(kznf1_pos - kznf2_pos) / 4000000
            return score
        else:
            return 0
    else:
        return -2

def defineKznfOrigin():
    '''
    Comprehensively consider the age of KZFPs,
    the distance between genes and the similarity information between KZFPs domains,
    and judge the source of KRAB and C2H2 domains of each KZFPs
    '''
    id_domainSeq_dic = np.load(r'\data\id_domainDone_dic.npy',allow_pickle=True)[()]
    KznfGeneAge_dic = np.load(r'\data\KznfGeneAge_dic.npy',allow_pickle=True)[()]
    EvolutionAge_dic = KznfGeneAge_dic['Evolution']
    geneId_dic = KznfGeneAge_dic['geneId']
    Idgene_dic = {}
    for keys,value in geneId_dic.items():
        Idgene_dic[value] = keys
    kznfDisAnno_df = pd.read_csv(r'\data\kznfDistance.csv')
    KznfDomainOrigin_dic = np.load(r'\data\KznfDomainOriginFromAllC2h2_AllKznf_40_dic.npy',
                                   allow_pickle=True)[()]
    KznfDomainOriginResult_dic = {}
    for kznf in KznfDomainOrigin_dic.keys():
        print(kznf)
        singleDomainResult_dic = {}
        singleKznfDomainAlign_dic = KznfDomainOrigin_dic[kznf]
        for domain in singleKznfDomainAlign_dic.keys():
            domainAlignInfro_list = singleKznfDomainAlign_dic[domain]
            domainScore_dic = {}
            for align in domainAlignInfro_list:
                Similirity_score = align[0]
                geneName = Idgene_dic[align[1].split('_')[0]]
                age = EvolutionAge_dic[geneName]
                aimKznf_age = EvolutionAge_dic[kznf]
                if age < aimKznf_age:
                    continue
                if age == aimKznf_age:
                    age = 0
                if age > aimKznf_age :
                    age = 2
                distence = getKznfDistanceScore(kznfDisAnno_df,kznf,geneName)
                combineScore = age + Similirity_score + distence
                domainScore_dic[combineScore] = [align[1],Similirity_score,distence]
            if domainScore_dic == {}:
                singleDomainResult_dic[domain] = {'id': geneId_dic[kznf],
                                                  'Similarity': 100,
                                                  'idEvoTime': EvolutionAge_dic[kznf],
                                                  'idDistence':0}
            else:
                maxCombineScore = max(list(domainScore_dic.keys()))
                idReuslt = domainScore_dic[maxCombineScore][0]
                singleDomainResult_dic[domain] = {'id':idReuslt,
                                                  'Similarity':domainScore_dic[maxCombineScore][1],
                                                  'idEvoTime':EvolutionAge_dic[Idgene_dic[idReuslt.split('_')[0]]],
                                                  'idDistence':domainScore_dic[maxCombineScore][2]}
        KznfDomainOriginResult_dic[kznf] = singleDomainResult_dic
    np.save(r'\data\KznfDomainOriginResult_All_dic.npy',KznfDomainOriginResult_dic)
    KznfDomainOriginResult_dic = np.load(r'\data\KznfDomainOriginResult_dic.npy',allow_pickle=True)[()]


