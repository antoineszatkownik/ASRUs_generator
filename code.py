import numpy as np
import pandas as pd
import csv
from copy import deepcopy
import math
import os
import sys
import itertools

sys.setrecursionlimit(100000)

eventsDup_path = 'C:/Users/Antoine/Desktop/BIM/projetS2/TranscriptAnnotation-master-duplications-curated_data/TranscriptAnnotation-master-duplications-curated_data/duplications/curated_data/eventsDupCons.txt'
names_path='C:/Users/Antoine/Desktop/BIM/projetS2/dupRaw'
allpath_path='C:/Users/Antoine/Desktop/BIM/projetS2/allPaths.txt'

allPaths=[]  # une ligne = gene + chemin can
with open(allpath_path,'r') as f:
    for lines in f.readlines():
        allPaths.append(lines.rstrip('\n').split(' '))

def get_immediate_subdirectories(a_dir): #while inside directory get subdirectories names into list
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def eventsDup_df(gene,path=eventsDup_path): #fonction qui return la portion dans eventsDup du gene targeté
    
    df_events=pd.read_csv(path)    #csv to pandas.dataframe
    df_events_gene=df_events.loc[df_events['gene']==gene]
    df_events_gene_np=df_events_gene.to_numpy()
    
    return df_events_gene_np

def ases_df(gene):  #fait le dataframe de l'ases du gene targete où on isole les colonnes d'intérets

    ases_gene_path= 'C:/Users/Antoine/Desktop/BIM/projetS2/aseDup/{}/ases_table.csv'.format(gene)
    df_ases_gene=pd.read_csv(ases_gene_path)[['CanonicalPath','AlternativePath','ASE','MutualExclusiveCanonical','MutualExclusiveAlternative']]
    df_ases_gene_np=df_ases_gene.to_numpy()

    return df_ases_gene_np

def sexSize(gene,sex):  #info sur le sex, sa taille est dans la colonne 3

    sexSizeEvents_path = 'C:/Users/Antoine/Desktop/BIM/projetS2/TranscriptAnnotation-master-duplications-curated_data/TranscriptAnnotation-master-duplications-curated_data/duplications/curated_data/sexSizeEvents.txt'
    df_size=pd.read_csv(sexSizeEvents_path)
    df_size_gene=df_size.loc[df_size['gene']==gene].loc[df_size['sexon']==sex].to_numpy()
    return df_size_gene

 
#fonction qui regroupe ensemble par transitivite les elements dune liste de listes
def transitive_closure(s):
    if len(s) < 2: return s
 
    r, b = [s[0]], transitive_closure(s[1:])
    for x in b:
        if r[0].intersection(x): r[0].update(x)
        else: r.append(x)
    return r

### A PARTIR D'ICI ON S'OCCUPE D'ETENDRE

def namestr(obj, namespace=globals()):  #get name of var
    return [name for name in namespace if namespace[name] is obj]

def event_to_pairs(gene): #{ b_i : [(paire,role),...,ligne dans ases], ... }

    dict={}

    eventsDup_df_np=eventsDup_df(gene)

    for rank in np.unique(eventsDup_df_np[:,3]):
        for row in np.where(eventsDup_df_np[:,3]==rank)[0]:
            if eventsDup_df_np[row][-1]!='No':
                pair,role=(eventsDup_df_np[row][1],eventsDup_df_np[row][2]),eventsDup_df_np[row][-1]
                if rank in dict:
                    dict[rank].append((pair,role))
                else:
                    dict[rank]=[(pair,role)]
        if rank in dict:
            dict[rank].append(ases_df(gene)[rank-1])

    return dict

def pair_to_events(gene): #{ (paire) : [b_1,...,b_p], ...}

    dict={}

    for row in eventsDup_df(gene):
        if row[-1]!='No':
            if tuple(row[1:3]) in dict:
                dict[tuple(row[1:3])].append(row[3])
            else:
                dict[tuple(row[1:3])]=[row[3]]

    return dict

def exons(gene,A,B):  #extrait l'info des hhr
    if A!=B:   
        try:
            aln='C:/Users/Antoine/Desktop/BIM/projetS2/dupRaw/{}/{}.{}.hhr'.format(gene,A,B)
        except FileNotFoundError:
            aln='C:/Users/Antoine/Desktop/BIM/projetS2/dupRaw/{}/{}.{}.hhr'.format(gene,B,A)

        info=[]
        f=open(aln,'r')
        for lines in f.readlines()[8:]:
            if lines.strip() != '':
                info.append(lines.split())
    else:
        return
    
    return np.array(info)

def extension_marge_pair(gene,A,B): #A est dans le can et B dans le alt

    aln=exons(gene,A,B)
    flag=False
    if '-' not in aln[1][-3]:  #gere des cas de bugs d'ecriture dans les hhr
        colA=aln[1][-2].split('-')
        flag=True
    else:
        colA=aln[1][-3].split('-')
    if flag==True and aln[1][-1].rsplit('(')[0].split('-')[0]!='':
        colB=aln[1][-1].rsplit('(')[0].split('-')
    else:
        colB=aln[1][-2].split('-')
    taille_A=float(aln[6][-1].strip('()'))
    try:
        taille_B=float(aln[1][-1].strip('()'))
    except ValueError:
        taille_B=float(aln[9][-1].strip('()'))

    coverage_A=np.abs(float(colA[1])-float(colA[0])+1)*100/taille_A
    coverage_B=np.abs(float(colB[1])-float(colB[0])+1)*100/taille_B

    if coverage_A<100 or coverage_B<100:
        if float(colB[0])==1:
            marge_B_Nter=0
            if float(colA[0])==1:
                marge_A_Nter=0
            if float(colA[0])>1:
                marge_A_Nter=float(colA[0])-1
        if float(colB[0])>1:
            marge_B_Nter=float(colB[0])-1
            if float(colA[0])==1:
                marge_A_Nter=0
            if float(colA[0])>1:
                marge_A_Nter=float(colA[0])-1
        marge_A_Cter=taille_A-float(colA[1])
        marge_B_Cter=taille_B-float(colB[1])
    else:
        marge_B_Nter, marge_B_Cter, marge_A_Nter, marge_A_Cter=0,0,0,0
        
        #si marge_B est positive alors on cest A quon etend sinon cest B

    return [marge_B_Nter, marge_B_Cter, marge_A_Nter, marge_A_Cter]
    
def concatenation(gene,A,C,direction):
    #direction == True means on etend A en Cter si cest False on l'etend en Nter
    
    path_A='C:/Users/Antoine/Desktop/BIM/projetS2/berenice/{}/msa_s_exon_{}.txt'.format(gene,A)
    path_C='C:/Users/Antoine/Desktop/BIM/projetS2/berenice/{}/msa_s_exon_{}.txt'.format(gene,C)

    with open(path_A,'r') as f:
        info_A = [line.rstrip('\n') for line in f]
    with open(path_C,'r') as f:  
        info_C = [line.rstrip('\n') for line in f]
    
    concat=[]
    
    if list(A)[0]==list(C)[0]:  #le nombre despeces et les especes sont les meme
        for doublet in list(zip(info_A,info_C)):
            if doublet[0][0]=='>':
                concat.append(doublet[0])
            else:
                if direction==True:
                    concat.append(doublet[0]+doublet[1])
                else:
                    concat.append(doublet[1]+doublet[0])
    #else  ajouter des lignes de gap


    f=open('msa_s_exon_{}.{}.txt'.format(A,C),'a')
    for line in concat:
        f.write(line+'\n')
    f.close()
            
    return np.array(concat)
    
#print(concatenation('TPM1','1_2','1_3',True))

def nettoyage(liste):  #enleve les doublons dans une unite
    singleton=[]
    for i in liste.copy():
        flag=False
        if '/' in i:
            continue
        else:
            for j in liste:
                if i==j:
                    continue
                else:
                    if '/' not in j:
                        continue
                    else:
                        if i in j.split('/'):
                            liste.remove(i)
                            flag=True
        if flag!=True:
            singleton.append(i)
            liste.remove(i)

    return set(liste).union(singleton)

def codage(eve2pair,C,d_X):   

    isgood=True
    for evebis in eve2pair:
        status_ext=(C in eve2pair[evebis][-1][0].split('/'),C in eve2pair[evebis][-1][1].split('/'),int(C in eve2pair[evebis][-1][0].split('/'))+int(C in eve2pair[evebis][-1][1].split('/')))
        print('-----------------------------',status_ext,d_X[evebis],evebis)
        if d_X[evebis][2]==2 and status_ext[2]==1:
            isgood=False
            break
        if d_X[evebis][2]==0 and status_ext[2]==1:
            isgood=False
            break
        if d_X[evebis][2]==1 and status_ext!=d_X[evebis]:
            isgood=False
            break

    return isgood

def getdata(gene,Y,Z): #recupere les infos dans les hhr et gere des cas de bugs d'ecriture dans les hhr
    try:
        aln=exons(gene,Y,Z) # on regarde l'alignement du candidat avec l'exon en face de celui qu'on etend
        if len(aln[1][-3].split('-'))==1:
            colY=aln[1][-2].split('-')
            if '(' in aln[1][-1]:
                colZ=aln[1][-1].split('(')[0].split('-')
            else:
                colZ=aln[1][-1].split('-')
        else:
            colY=aln[1][-3].split('-')
            colZ=aln[1][-2].split('-')
        taille_Y=float(aln[5][-1].strip('()'))
        try:
            taille_Z=float(aln[1][-1].strip('()'))
        except ValueError:
            taille_Z=float(aln[1][-1].rsplit('(')[1].strip(')'))
    except FileNotFoundError:
        aln=exons(gene,Z,Y)
        if len(aln[1][-3].split('-'))==1:
            colZ=aln[1][-2].split('-')
            if '(' in aln[1][-1]:
                colY=aln[1][-1].split('(')[0].split('-')
            else:
                colY=aln[1][-1].split('-')
        else:
            colZ=aln[1][-3].split('-')
            colY=aln[1][-2].split('-')
        taille_Z=float(aln[5][-1].strip('()'))
        try:
            taille_Y=float(aln[1][-1].strip('()'))
        except ValueError:
            taille_Y=float(aln[1][-1].rsplit('(')[1].strip(')'))

    return colY,colZ,taille_Y,taille_Z

def most_inner_loop(gene,pair,pair2eve,pair2evecopy,eve2pair,allPaths,extension_marge,marge,d_A,d_B,A,B):
    for eve in pair2evecopy[(A,B)]:
        print(eve,pair)
        if extension_marge.index(marge)<=1:  #path_bool me dit si je regarde le chemin can ou alt
            path_bool=0  #chemin can
            X=A  #on etend A
            d_X=d_A
            if extension_marge.index(marge)==0:
                direction='Nter'
            else:
                direction='Cter'
        else:
            path_bool=1  #chemin alt
            X=B  #on etend B
            d_X=d_B
            if extension_marge.index(marge)==2:
                direction='Nter'
            else:
                direction='Cter'

        path=eve2pair[eve][-1][path_bool].split('/') #je recupere le chemin can ou alt
        if X not in path:
            path=allPaths[np.where(np.array(allPaths)[:,0]==gene)[0][0]][1].split('/')
            if X not in path:
                continue

        if len(path)>1 and ((direction=='Cter' and path.index(X)<len(path)-1) or (direction=='Nter' and path.index(X)>0)):
            if direction=='Cter':
                C=path[path.index(X)+1]  # en Cter le candidat est le noeud pointe par X
            else:
                C=path[path.index(X)-1]  #en Nter le candidat est le noeud qui pointe X

            if C=='start' or C=='stop' or (C not in eve2pair[eve][-1][path_bool].split('/')):
                continue
            else:
                if C.split('_')[0]=='0':
                    continue

                else:
                    if sexSize(gene,C)[0][3]<=5: #si le candidat est de taille inf a 5 je l'ajoute
                        print('extension inf a 5',eve,pair,X,C,direction)
                        flag=False
                        isgood=codage(eve2pair,C,d_X)

                        if isgood:
                            if direction=='Cter' and X==A:
                                update=(X+'/'+C,B)
                            elif direction=='Cter' and X==B:
                                update=(A,X+'/'+C)
                            elif direction=='Nter' and X==A:
                                update=(C+'/'+X,B)
                            else:
                                update=(A,C+'/'+X)
                            pair2eve[update]=pair2eve.pop(pair,pair2evecopy[pair])
                        else:
                            continue

                    else:
                        aln_bis=exons(gene,A,B) #on regarde l'alignement de la pair
                        taille_A=float(aln_bis[5][-1].strip('()'))
                        print(A,B)
                        if len(aln_bis[1][-3].split('-'))==1:
                            colA_bis=aln_bis[1][-2].split('-')
                        else:
                            colA_bis=aln_bis[1][-3].split('-')
                        if '' in aln_bis[1][-1].split('('):
                            colB_bis=aln_bis[1][-2].split('-')
                        else:
                            colB_bis=aln_bis[1][-1].split('(')[0].split('-')

                        if X==A:
                            if B!=C:
                                get_data=getdata(gene,C,B)
                                colB=get_data[1]
                                if (float(colB_bis[0])<=float(colB[0])<=float(colB_bis[1]) or float(colB_bis[0])<=float(colB[1])<=float(colB_bis[1])) or (float(colB[0])<=float(colB_bis[0])<=float(colB[1]) or float(colB[0])<=float(colB_bis[1])<=float(colB[1])):
                                    continue #on verifie que lalignement entre A et B et l'alignement correspondant impliquant l'extension ne s'overlap pas
                                else:
                                    isgood=codage(eve2pair,C,d_X)

                                    if isgood:
                                        if direction=='Cter':
                                            if float(colB_bis[1])<float(colB[0]): #on verifie que les bornes de lalignements sont coherentes
                                                update=(X+'/'+C,B)
                                            else:
                                                continue

                                        elif direction=='Nter':
                                            if float(colB[1])<float(colB_bis[0]):
                                                update=(C+'/'+X,B)
                                            else:
                                                continue
                                    else:
                                        continue
                                    pair2eve[update]=pair2eve.pop(pair,pair2evecopy[pair])

                            else:
                                continue
                        else:
                            if A!=C:
                                get_data=getdata(gene,A,C)
                                colA=get_data[0]
                                if (float(colA_bis[0])<=float(colA[0])<=float(colA_bis[1]) or float(colA_bis[0])<=float(colA[1])<=float(colA_bis[1])) or (float(colA[0])<=float(colA_bis[0])<=float(colA[1]) or float(colA[0])<=float(colA_bis[1])<=float(colA[1])):
                                    continue
                                else:
                                    isgood=codage(eve2pair,C,d_X)

                                    if isgood:
                                        if direction=='Cter':
                                            if float(colA_bis[1])<float(colA[0]):
                                                update=(X+'/'+C,A)  #ITOU WARNING
                                            else:
                                                continue

                                        if direction=='Nter':
                                            if float(colA[1])<float(colA_bis[0]):
                                                update=(C+'/'+X,A)    #ATTENTION ICI AVANT A CT B
                                            else:
                                                continue

                                    else:
                                        continue
                                    pair2eve[update]=pair2eve.pop(pair,pair2evecopy[pair])

                            else:
                                continue

def extensionbis(gene):

    eve2pair=event_to_pairs(gene)
    pair2eve=pair_to_events(gene)
    print(pair2eve)

    pair2evecopy=deepcopy(pair2eve)

    for elem in list(pair2eve):  #on retire les paires avec identite < 1%
        if float(exons(gene,elem[0],elem[1])[4][4].split('=')[1].strip('%'))<=1:
            del pair2eve[elem]
        else:
            continue

    #on applique la transitivite sur les paires
    temp=transitive_closure([set(elem) for elem in list(pair2eve)])

    #on recupere les evenements pour lecriture du csv
    events_unit=[]
    for unit in temp:
        events_for1RU_tmp=[]
        for pair in pair2eve:
            if pair[0] in unit and pair[1] in unit:
                try:
                    events_for1RU_tmp.append(pair2eve[(pair[0],pair[1])])
                except KeyError:
                    events_for1RU_tmp.append(pair2eve[(pair[1],pair[0])])
            else:
                continue     
        events_unit.append(np.unique(list(itertools.chain.from_iterable(events_for1RU_tmp))))

    for pair in pair2eve.copy():  #je loop sur toutes les paires 
        A,B=pair[0],pair[1]
        d_A={}
        d_B={}
        if float(exons(gene,A,B)[4][4].split('=')[1].strip('%'))<=1: #si le pourcentage didentite est inf a 1 alors faux positif
            del pair2eve[pair]
        else:
            extension_marge=extension_marge_pair(gene,A,B)  #je calcul si je peux etendre 
            print(extension_marge, pair)

            if extension_marge[0]>0 or extension_marge[1]>0:
                for eve in eve2pair:
                    d_A[eve]=(A in eve2pair[eve][-1][0].split('/'),A in eve2pair[eve][-1][1].split('/'),int(A in eve2pair[eve][-1][0].split('/'))+int(A in eve2pair[eve][-1][1].split('/')))
            if extension_marge[2]>0 or extension_marge[3]>0:
                for eve in eve2pair:
                    d_B[eve]=(B in eve2pair[eve][-1][0].split('/'),B in eve2pair[eve][-1][1].split('/'),int(B in eve2pair[eve][-1][0].split('/'))+int(B in eve2pair[eve][-1][1].split('/')))
            print('-----------',d_A,d_B)
            for marge in extension_marge:
                if marge<=0:
                    continue
                else:
                    most_inner_loop(gene,pair,pair2eve,pair2evecopy,eve2pair,allPaths,extension_marge,marge,d_A,d_B,A,B)

    for unite in temp:  #gere les cas d'overlap entre instances 
        flag=False
        for pair in list(pair2eve):
            if '/' in pair[0]:
                if pair[0].split('/')[0] in unite and pair[1] in unite:
                    unite.remove(pair[0].split('/')[0])
                    unite.add(pair[0])
                elif pair[0].split('/')[1] in unite and pair[1] in unite:
                    unite.remove(pair[0].split('/')[1])
                    unite.add(pair[0])
                else:
                    continue
            elif '/' in pair[1]:
                if pair[1].split('/')[0] in unite and pair[0] in unite:
                    unite.remove(pair[1].split('/')[0])
                    unite.add(pair[1])
                elif pair[1].split('/')[1] in unite and pair[0] in unite:
                    unite.remove(pair[1].split('/')[1])
                    unite.add(pair[1])
                else:
                    continue

    for instance in temp:    #itou
        for seed in instance.copy():
            if '/' in seed:
                if seed.split('/')[0] in instance:
                    instance.remove(seed.split('/')[0])
                elif seed.split('/')[1] in instance:
                    instance.remove(seed.split('/')[1])
            else:
                continue

    tempbis=[]
    for inst in temp:
        tempbis.append(nettoyage(list(inst)))


    for unit in tempbis:   #sert a faire des instances a 3 s-exons
        for couple in list(itertools.combinations(list(unit),2)):
            if couple[0] in unit and couple[1] in unit:
                if '/' in couple[0] and '/' in couple[1]:
                    if couple[0].split('/')[0]==couple[1].split('/')[1]:
                        unit.remove(couple[0])
                        unit.remove(couple[1])
                        unit.add(couple[1].split('/')[0]+'/'+couple[0])
                    if couple[0].split('/')[1]==couple[1].split('/')[0]:
                        unit.remove(couple[0])
                        unit.remove(couple[1])
                        unit.add(couple[0]+'/'+couple[1].split('/')[1])
                else:
                    continue
            else:
                continue

    print('----------------------------------------------')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('----------------------------------------------')

    return tempbis,pair2eve,events_unit

#extensionbisdata=extensionbis('ANK2')
#print(extensionbisdata)


earth=get_immediate_subdirectories(names_path)
#POUR ECRIRE UNE TABLE AVEC 1 LIGNE = 1 UNITE REPETEE
def writecsv_ASRU():

    with open('datatest__writeNEWW.csv', 'w',newline="") as new_file:
        csv_writer=csv.writer(new_file, delimiter=',')

        csv_writer.writerow(['gene','uniteRepetee','#instances','max','min','moy','median','ecartType','evenements'])

        for gene in earth:
            print(gene)
            print('--------------------------------------------')
            print('--------------------------------------------')
            
            data=extensionbis(gene)
            RU_genedata=data[0]
            pair2eve_updated=data[1]
            events_unit=data[2]

            for unit in RU_genedata:
                gene_csv_line=[]
                gene_csv_line.append(gene)   #nom du gene
                gene_csv_line.append(unit)
                gene_csv_line.append(len(unit))
                instancedatatmp=[]

                #retrouver les paires qui ont permis de construire la RU, recuperer les evenements des paires dune RU
                # events_for1RU_tmp=[]
                # for pair in list(itertools.combinations(list(unit), 2)):
                #     if (pair[0],pair[1]) in list(pair2eve_updated) or (pair[1],pair[0]) in list(pair2eve_updated):
                #         try:
                #             events_for1RU_tmp.append(pair2eve_updated[(pair[0],pair[1])])
                #         except KeyError:
                #             events_for1RU_tmp.append(pair2eve_updated[(pair[1],pair[0])])
                #     else:
                #         continue

                
                # events_for1RU=np.unique(list(itertools.chain.from_iterable(events_for1RU_tmp)))

                for inst in unit:
                    if '/' not in inst:
                        instancedatatmp.append(sexSize(gene,inst)[0][3])
                    elif len(inst.split('/'))==2:
                        instancedatatmp.append(sexSize(gene,inst.split('/')[0])[0][3]+sexSize(gene,inst.split('/')[1])[0][3])
                    elif len(inst.split('/'))==3:
                        instancedatatmp.append(sexSize(gene,inst.split('/')[0])[0][3]+sexSize(gene,inst.split('/')[1])[0][3]+sexSize(gene,inst.split('/')[2])[0][3])
                gene_csv_line.append(max(instancedatatmp))
                gene_csv_line.append(min(instancedatatmp))
                gene_csv_line.append(np.mean(instancedatatmp))
                gene_csv_line.append(np.median(instancedatatmp))
                gene_csv_line.append(np.std(instancedatatmp))
                gene_csv_line.append(events_unit[RU_genedata.index(unit)])
                csv_writer.writerow(gene_csv_line)

#writecsv_ASRU()

# POUR ECRIRE UNE TABLE AVEC 1LIGNE=1INSTANCE
def writecsv_instances():
    with open('datastats__instNEWW.csv', 'w',newline="") as new_file:
        csv_writer=csv.writer(new_file, delimiter=',')

        csv_writer.writerow(['instance','taille','nombre','UniteRepetee','gene'])

        for gene in earth:
            print(gene)
            print('--------------------------------------------')
            #print('les unites repetees avant extension {}'.format(instanciation('ANK2',eventsDup_df('ANK2'))))

            print('--------------------------------------------')
            
            data=extensionbis(gene)
            RU_genedata=data[0]
            pair2eve_updated=data[1]

            for unit in RU_genedata:
                #unit2pair=list(itertools.combinations(list(unit), 2))
                dejavu=[]
                for inst in unit:
                    if inst in dejavu:
                        continue
                    else:
                        gene_csv_line=[]
                        gene_csv_line.append(inst)
                        if '/' not in inst:
                            gene_csv_line.append(sexSize(gene,inst)[0][3])
                            gene_csv_line.append(1)
                        elif len(inst.split('/'))==2:
                            gene_csv_line.append(sexSize(gene,inst.split('/')[0])[0][3]+sexSize(gene,inst.split('/')[1])[0][3])
                            gene_csv_line.append(2)
                        else:
                            gene_csv_line.append(sexSize(gene,inst.split('/')[0])[0][3]+sexSize(gene,inst.split('/')[1])[0][3]+sexSize(gene,inst.split('/')[2])[0][3])
                            gene_csv_line.append(3) 

                        gene_csv_line.append(unit)
                        gene_csv_line.append(gene)
                        csv_writer.writerow(gene_csv_line)
                        dejavu.append(inst)

            else:
                continue

#writecsv_instances()