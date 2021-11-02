import numpy as np
import pandas as pd
import csv
from copy import deepcopy
import math
import os
import sys
import itertools
import networkx as nx
import matplotlib.pyplot as plt

sys.setrecursionlimit(100000)


eventsDup_path = '/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/eventsDupCons.txt'
names_path='/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/dupRaw'
allpath_path='/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/allPaths.txt'
proteomeDup_path='/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/proteome_duplication_pairs.csv'


allPaths=[]  # une ligne = gene + chemin can entier
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
    
    return df_events_gene_np,df_events_gene

def ases_df(gene):  #fait le dataframe de l'ases du gene targete où on isole les colonnes d'intérets (chemin can, alt, ...) voir 2eme variable

    ases_gene_path= '/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/aseDup/{}/ases_table.csv'.format(gene)
    df_ases_gene=pd.read_csv(ases_gene_path)[['CanonicalPath','AlternativePath','ASE','MutualExclusiveCanonical','MutualExclusiveAlternative']] #ICI
    df_ases_gene_np=df_ases_gene.to_numpy()

    return df_ases_gene_np

def sexSize(gene,sex):  #info sur le sex, sa taille est dans la colonne 3

    sexSizeEvents_path = '/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/sexSizeEvents.txt'
    df_size=pd.read_csv(sexSizeEvents_path)
    df_size_gene=df_size.loc[df_size['gene']==gene].loc[df_size['sexon']==sex].to_numpy()
    return df_size_gene


# def sim_graph(gene): #construction du graphe de similarité d'un gène, on enlève les paires avec id_align <= 10% et les paires de type 'NO'
#     proteomeDup_df=pd.read_csv(proteomeDup_path)
#     proteomeDup_gene=proteomeDup_df.loc[proteomeDup_df['Gene']==gene]
#     proteomeDup_gene=proteomeDup_gene.loc[proteomeDup_gene['Identity']>10].to_numpy()[:,0:2]

#     #si on veut éliminer les paires UNREL
#     #eventsDup_df_NO=eventsDup_df(gene)[1].loc[eventsDup_df(gene)[1]["typePair"] != 'No']
#     #eventsDup_df_NO_UNREL=list(eventsDup_df_NO.loc[eventsDup_df_NO["typePair"]!='UNREL'].to_numpy()[:,1:3])

#     #si on veut juste juste éliminer NO
#     eventsDup_df_NO=list(eventsDup_df(gene)[1].loc[eventsDup_df(gene)[1]["typePair"] != 'No'].to_numpy()[:,1:3])

#     return nx.from_edgelist(np.array([x for x in set(tuple(x) for x in eventsDup_df_NO) & set(tuple(x) for x in proteomeDup_gene)]), create_using=nx.Graph)

#nx.draw(sim_graph('AHNAK'),with_labels=True)
#plt.show()

def sim_graph(gene): #construction du graphe de similarité d'un gène, on enlève les paires avec id_align <= 10% et les paires de type 'NO'
    proteomeDup_df=pd.read_csv(proteomeDup_path)
    proteomeDup_gene=proteomeDup_df.loc[proteomeDup_df['Gene']==gene]
    proteomeDup_gene=proteomeDup_gene.loc[proteomeDup_gene['Identity']>10].to_numpy()[:,0:2]
    eventsDup_gene=eventsDup_df(gene)[1].loc[eventsDup_df(gene)[1]["typePair"] != 'No'].to_numpy()[:,1:3]
#     #si on veut éliminer les paires UNREL et remplacer dans le return eventsDup_gene par eventsDup_df_NO_UNREL
#     #eventsDup_df_NO=eventsDup_df(gene)[1].loc[eventsDup_df(gene)[1]["typePair"] != 'No']
#     #eventsDup_df_NO_UNREL=eventsDup_df_NO.loc[eventsDup_df_NO["typePair"]!='UNREL'].to_numpy()[:,1:3]
    return nx.from_edgelist([pair for pair in eventsDup_gene if pair in proteomeDup_gene], create_using=nx.Graph)

#ici les paires viennent directement du graphe de similarité, donc une modification sur le graphe doit changer
#ces 2 dicos

def event_to_pairs(gene):

    #{ b_i : [paire,...,ligne dans ases], ... } dicou où les clés sont les ranks des évènements
    #les valeurs sont les paires impliquées dans cet évènement et les chemins can et alt de l'evt

    g=sim_graph(gene)
        
    dict={}

    eventsDup_df_np=eventsDup_df(gene)[0]

    for rank in np.unique(eventsDup_df_np[:,3]): #on boucle sur les evts
        for row in np.where(eventsDup_df_np[:,3]==rank)[0]: #on boucle sur les paires impliquées dans cet evt
            if eventsDup_df_np[row][-1]!='No' and ((eventsDup_df_np[row][1],eventsDup_df_np[row][2]) in list(g.edges()) or (eventsDup_df_np[row][2],eventsDup_df_np[row][1]) in list(g.edges()) ):
                try:
                    pair_in_g=list(g.edges())[list(g.edges()).index((eventsDup_df_np[row][1],eventsDup_df_np[row][2]))]
                except ValueError:
                    pair_in_g=list(g.edges())[list(g.edges()).index((eventsDup_df_np[row][2],eventsDup_df_np[row][1]))]
                if rank in dict:
                    dict[rank].append(pair_in_g)
                else:
                    dict[rank]=[pair_in_g]
        if rank in dict:
            dict[rank].append(ases_df(gene)[rank-1])

    return dict,g

def pair_to_events(gene): 

    #{ (paire) : ([b_1,...,b_p],classe), ...} dico où les clés sont les paires
    #les valeurs sont les évènements dans lesquels cette paire est impliquée et la classe (MEX,ALT,REL,UNREL) de cette paire

    g=sim_graph(gene)
    
    dict={}
    
    eventsDup_df_np=eventsDup_df(gene)[0]

    for row in eventsDup_df_np: #on boucle sur les paires
        if row[-1]!='No' and ((row[1],row[2]) in list(g.edges()) or (row[2],row[1]) in list(g.edges())):
            role=row[-1]
            try:
                pair_in_g=list(g.edges())[list(g.edges()).index((row[1],row[2]))]
            except ValueError:
                pair_in_g=list(g.edges())[list(g.edges()).index((row[2],row[1]))]
            if tuple(row[1:3]) in dict:
                dict[pair_in_g][0][0].append(row[3])
            else:
                dict[pair_in_g]=[([row[3]],role)]
            

    return dict

def exons(gene,A,B):  #extrait l'info des hhr
    if A!=B:   
        info=[]
        try:
            aln='/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/dupRaw/{}/{}.{}.hhr'.format(gene,A,B)
            f=open(aln,'r')
        except FileNotFoundError:
            aln='/Users/antoineszatkownik/Documents/projetAS/TranscriptAnnotation/duplications/curated_data/dupRaw/{}/{}.{}.hhr'.format(gene,B,A)
            f=open(aln,'r')
        for lines in f.readlines()[8:]:
            if lines.strip() != '':
                info.append(lines.split())
    else:
        return
    
    return np.array(info)

def extension_marge_pair(gene,A,B): 

    #calcul pour un alignement entre A et B si je peux étendre A ou B en N ou C terminal

    #A est dans le can et B dans le alt

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

def codage(eve2pair,C,d_X): 

    #C est le candidat pour faire l'extension
    #voir pdf du rapport pour les détails du codage

    isgood=True
    for evebis in eve2pair: #on boucle sur les évènements pour vérifier si l'extension est compatible avec tout les évènements
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

def getdata(gene,Y,Z): 
    
    #recupere les infos dans les hhr (alignement entre le candidat et le sex en face de celui qu'on étend)  
    #et gere des cas de bugs d'ecriture dans les hhr

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

def most_inner_loop(gene,g,pair,pair2eve,pair2evecopy,eve2pair,allPaths,extension_marge,marge,d_A,d_B,A,B):

    #fonction principale qui fait les extensions et update les noeuds du graphe

    for eve in pair2evecopy[(A,B)][0][0]: #pour chaque paire on boucle sur ses évènements
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
        #X est celui qu'on étend
        path=eve2pair[eve][-1][path_bool].split('/') #je recupere le chemin can ou alt
        if X not in path: #si X n'est pas dans le chemin can ou alt, cad que c'est une ancre
            path=allPaths[np.where(np.array(allPaths)[:,0]==gene)[0][0]][1].split('/') #je recupere le chemin canonique global
            if X not in path:
                continue

        #on vérifie que la taille du chemin et qu'il y a bien des éléments avant ou après X
        #C est le candidat pour faire l'extension (qui sera de la forme C/X ou X/C)
        if len(path)>1 and ((direction=='Cter' and path.index(X)<len(path)-1) or (direction=='Nter' and path.index(X)>0)): 
            if direction=='Cter':
                C=path[path.index(X)+1]  # en Cter le candidat est le noeud pointe par X
            else:
                C=path[path.index(X)-1]  #en Nter le candidat est le noeud qui pointe X

            if C=='start' or C=='stop' or (C not in eve2pair[eve][-1][path_bool].split('/')): #check if C is not start or stop node of ESG
                continue
            else:
                if C.split('_')[0]=='0': #check if C is not a node '0_blabla' which are not real sex
                    continue

                else:
                    if sexSize(gene,C)[0][3]<=5: #si le candidat est de taille inf a 5 je l'ajoute
                        print('extension inf a 5',eve,pair,X,C,direction)
                        flag=False
                        isgood=codage(eve2pair,C,d_X)

                        if isgood:
                            #on update g.nodes()
                            if direction=='Cter' and (X==A or X==B):
                                update=X+'/'+C
                            elif direction=='Nter' and (X==A or X==B):
                                update=C+'/'+X
                            mapping={X:update}
                            print(mapping,'MAPPING')
                            g=nx.relabel_nodes(g, mapping)  #j'update le noeud seed X du graphe en le noeud seed-candidat
                            try:  # l'unite du type { X , C , ...} devient { X/C , C , ...} après extension, il reste à enlever le doublon C
                                g.remove_node(C)
                            except nx.exception.NetworkXError: # cas où C a déjà été enlevé
                                pass
                        else:
                            continue

                    else: #si le noeud est de taille >= 5 je dois regarder calculer les hhr de X et celui en face de X et du candidat avec celui en face de X
                        aln_bis=exons(gene,A,B) #on regarde l'alignement de la paire (seed, sex en face de seed)
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

                        if X==A: #dans la paire (A,B) la seed X est A
                            if B!=C:
                                get_data=getdata(gene,C,B) #on regarde l'alignement de la paire (candidat, sex en face de seed)
                                colB=get_data[1]
                                if (float(colB_bis[0])<=float(colB[0])<=float(colB_bis[1]) or float(colB_bis[0])<=float(colB[1])<=float(colB_bis[1])) or (float(colB[0])<=float(colB_bis[0])<=float(colB[1]) or float(colB[0])<=float(colB_bis[1])<=float(colB[1])):
                                    continue #on verifie que lalignement entre A et B et l'alignement correspondant impliquant l'extension ne s'overlap pas
                                else:
                                    isgood=codage(eve2pair,C,d_X)

                                    if isgood:
                                        if direction=='Cter':
                                            if float(colB_bis[1])<float(colB[0]): #on verifie que les bornes de lalignement sont coherentes
                                                update=X+'/'+C
                                            else:
                                                continue

                                        elif direction=='Nter':
                                            if float(colB[1])<float(colB_bis[0]):
                                                update=C+'/'+X
                                            else:
                                                continue
                                    else:
                                        continue
                                    mapping={X:update}
                                    print(mapping,'MAPPING')
                                    g=nx.relabel_nodes(g, mapping)
                                    try:
                                        g.remove_node(C)
                                    except nx.exception.NetworkXError:
                                        pass

                            else:
                                continue
                        else: #dans la paire (A,B) la seed X est B
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
                                                update=X+'/'+C  #ITOU WARNING
                                            else:
                                                continue

                                        if direction=='Nter':
                                            if float(colA[1])<float(colA_bis[0]):
                                                update=C+'/'+X    #ATTENTION ICI AVANT A CT B
                                            else:
                                                continue

                                    else:
                                        continue
                                    mapping={X:update}
                                    print(mapping,'MAPPING')
                                    g=nx.relabel_nodes(g, mapping)
                                    try:
                                        g.remove_node(C)
                                    except nx.exception.NetworkXError:
                                        pass

                            else:
                                continue
    return g #retourne le graphe avec les noeuds updatés

def extensionbis(gene):

    eve2pair=event_to_pairs(gene)[0]  #dico des evenements aux paires
    pair2eve=pair_to_events(gene)     #dico des paires aux évenements
    g=event_to_pairs(gene)[1]         #creation du squelette du sim graph (squelette car aucun extension n'est faite)
    print(pair2eve)

    #on fait une copie du dico car on ne peut pas boucler sur une structure qu'on modifie en même temps
    #mais cette variable est inutile dans cette version du code (pas sûr)
    pair2evecopy=deepcopy(pair2eve)

    #dans cette nouvelle version on devrait boucler sur les arêtes du graphe ??
    #plutot que le dico pair2eve. Avant pair2eve était la structure updatée maintenant c'est les nodes du graphes
    for elem in list(pair2eve):  #on retire les paires avec identite < 1%
        #if elem not in list(g.nodes())
        if float(exons(gene,elem[0],elem[1])[4][4].split('=')[1].strip('%'))<=1:
            del pair2eve[elem]
        else:
            continue

    #on applique la transitivite sur les paires

    #on recupere les evenements pour lecriture du csv
    # events_unit=[]
    # for unit in temp:
    #     events_for1RU_tmp=[]
    #     for pair in pair2eve:
    #         if pair[0] in unit and pair[1] in unit:
    #             try:
    #                 events_for1RU_tmp.append(pair2eve[(pair[0],pair[1])])
    #             except KeyError:
    #                 events_for1RU_tmp.append(pair2eve[(pair[1],pair[0])])
    #         else:
    #             continue     
    #     events_unit.append(np.unique(list(itertools.chain.from_iterable(events_for1RU_tmp))))

    #pourquoi refaire une copie ?
    for pair in pair2eve.copy():  #je loop sur toutes les paires 
        A,B=pair[0],pair[1]
        d_A={}
        d_B={}
        if float(exons(gene,A,B)[4][4].split('=')[1].strip('%'))<=1: #si le pourcentage didentite est inf a 1 alors faux positif
            del pair2eve[pair]
        else:
            extension_marge=extension_marge_pair(gene,A,B)  #je calcul si je peux etendre 
            print(extension_marge, pair)

            if extension_marge[0]>0 or extension_marge[1]>0: #Si j'étend A
                for eve in eve2pair:
                    d_A[eve]=(A in eve2pair[eve][-1][0].split('/'),A in eve2pair[eve][-1][1].split('/'),int(A in eve2pair[eve][-1][0].split('/'))+int(A in eve2pair[eve][-1][1].split('/')))
            if extension_marge[2]>0 or extension_marge[3]>0: #si j'étend B
                for eve in eve2pair:
                    d_B[eve]=(B in eve2pair[eve][-1][0].split('/'),B in eve2pair[eve][-1][1].split('/'),int(B in eve2pair[eve][-1][0].split('/'))+int(B in eve2pair[eve][-1][1].split('/')))
            print('-----------',d_A,d_B)
            for marge in extension_marge: #je considère l'extension pour toutes les marges non nulles
                if marge<=0:
                    continue
                else:
                    g=most_inner_loop(gene,g,pair,pair2eve,pair2evecopy,eve2pair,allPaths,extension_marge,marge,d_A,d_B,A,B) #j'update le squelette du sim graphe


    nx.draw(g, with_labels=True)
    plt.show()
    return list(nx.connected_components(g))

print(extensionbis('AHNAK'))

