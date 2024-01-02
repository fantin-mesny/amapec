import os
import sys
import argparse
import pandas as pd
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectFdr, chi2
from sklearn.calibration import CalibratedClassifierCV
import joblib


def get_params(argv):
    parser = argparse.ArgumentParser(description='Train ML model based on sequence properties, structure properties and k-mer content')
    parser.add_argument('-o', '--o', help="Output directory", required=True)
    parser.add_argument('-seq', '--seq', help="CSV file with sequence properties", required=True)
    parser.add_argument('-struc', '--struc', help="CSV file with structure properties", required=True)
    parser.add_argument('-kmer', '--kmer', help="CSV file with Kmer content", required=True)
#    parser.add_argument('-k', '--k', help="K size to use", default='3,4,5,6')
    parser.add_argument('-posMeta', '--posMeta', help="Positive metadata csv file", required=True)
    parser.add_argument('-negMeta', '--negMeta', help="Negative metadata csv file", required=True)
    parser.add_argument('-posPrefix','--posPrefix', help="Comma-separated ID prefixes (one letter) identifying positive training set",required=True)
    parser.add_argument('-negPrefix','--negPrefix', help="Comma-separated ID prefixes (one letter) identifying negative training set",required=True)
    parser.add_argument('--svmScale', action=argparse.BooleanOptionalAction)
    parser.add_argument('--excludeKmers', action=argparse.BooleanOptionalAction)
    parser.add_argument('--weighted', action=argparse.BooleanOptionalAction)
    a = parser.parse_args()
    return a

def getModelQC(preds):
    R2=preds['Correct prediction'].sum()/len(preds)
    Precision=preds['True positive'].sum()/(preds['True positive'].sum()+preds['False positive'].sum())
    Recall=preds['True positive'].sum()/(preds['True positive'].sum()+preds['False negative'].sum()) # aka Sensitivity aka True positive rate
    Specificity=preds['True negative'].sum()/(preds['False positive'].sum()+preds['True negative'].sum())
    BalancedR2=(Recall+Specificity)/2
    F1=(2*Precision*Recall)/(Precision+Recall)
    return {'R2':R2,'Precision':Precision,'Recall':Recall,'Specificity':Specificity,'Balanced R2':BalancedR2,'F1':F1}


if __name__ == '__main__':

    a = get_params(sys.argv[1:])
    
    ####### Parse inputs, combine properties into a single dataframe ######
    seq=pd.read_csv(a.seq).set_index('Name').drop(columns='Unnamed: 0').sort_index()
    struc=pd.read_csv(a.struc).set_index('Unnamed: 0').sort_index()
    
    if a.excludeKmers:
        allProp=pd.concat([seq,struc],axis=1).dropna()
    else:
        kmer=pd.read_csv(a.kmer).set_index('Unnamed: 0').sort_index()
#        kmer=kmer[[c for c in kmer.columns if len(c) in [int(val) for val in a.k.split(',')]]]
        kmer[kmer>1]=1
        allProp=pd.concat([seq,struc,kmer],axis=1).dropna()
        

    ####### Classify into negative and positive training set from prefix #######
    posPrefixes=str(a.posPrefix).split(',')
    negPrefixes=str(a.negPrefix).split(',')
    for ind in allProp.index:
        if ind[0] in posPrefixes:
            allProp.loc[ind,'Group']=1
        elif ind[0] in negPrefixes:
            allProp.loc[ind,'Group']=0
        else:
            print('ERROR: unrecognized prefix',ind[0])

    print('Number of proteins in training dataset:', len(allProp))
    print('        - in positive training dataset:', len(allProp[allProp['Group']==1]))
    print('        - in negative training dataset:', len(allProp[allProp['Group']==0]))
    
   
    ############ Filtering k-mers ############
    if not a.excludeKmers:
        print('Total number of K-mers counted:',len(kmer.columns))
        chi2_filter = SelectFdr(chi2,alpha=0.05)  #only keep enriched/depleted ones in positive vs. negative set according to chi2 test
        chi2_filter.fit(allProp[[c for c in kmer.columns]],allProp['Group'])
        kept=list(chi2_filter.get_feature_names_out())
        print('Number of K-mers kept after chi2:',len(kept))
        allProp=allProp[list(seq.columns)+list(struc.columns)+kept+['Group']]
   
    allProp.to_csv(a.o+'/trainingDataset.csv')
    print('Number of features in model:',len(allProp.columns)-1)
    
    
    ################ Standard scaling of the data ##################
    scaler = StandardScaler()
    allProp_scaled = pd.DataFrame(scaler.fit_transform(allProp.drop(columns='Group')))
    allProp_scaled.index=allProp.index
    allProp_scaled.columns=allProp.columns[:-1]
    allProp_scaled=allProp_scaled.merge(allProp[['Group']], left_index=True, right_index=True)
    allProp_scaled.to_csv(a.o+'/trainingDataset_scaled.csv')
    
    # export scale for later classification
    scaling=pd.DataFrame(index=allProp.columns,columns=['Standard deviation','Mean'])
    scaling['Standard deviation']=allProp.std(axis=0)
    scaling['Mean']=allProp.mean(axis=0)

    if a.svmScale:
        data=allProp_scaled
        scaling.to_csv(a.o+'/classifier.sav.scale')
    else:
        data=allProp
    ################### Perform leave-one-out CV ###################
    if a.weighted:
        weights={0: 1, 1: 2}
    else:
        weights={0:1,1:1}
    
    preds=[]
    for ind in data.index:
        dataOO=data.drop(index=[ind])
        clf = svm.SVC(kernel='linear',class_weight=weights)
        clf.fit(dataOO.drop(columns=['Group']), dataOO['Group'])

        pred=pd.DataFrame(columns=clf.classes_, index=[ind])
        predi=clf.predict(data[data.index==ind].drop(columns=['Group']))
        pred.loc[ind,'Prediction']=predi
        if predi==1 and ind[0] in posPrefixes:
        	pred.loc[ind,'True positive']=True
        	pred.loc[ind,'Correct prediction']=True
        elif predi==0 and ind[0] in negPrefixes:
        	pred.loc[ind,'True negative']=True
        	pred.loc[ind,'Correct prediction']=True
        elif predi==1 and ind[0] in negPrefixes:
        	pred.loc[ind,'False positive']=True 
        	pred.loc[ind,'Correct prediction']=False
        elif predi==0 and ind[0] in posPrefixes:
        	pred.loc[ind,'False negative']=True
        	pred.loc[ind,'Correct prediction']=False       	
        preds.append(pred)
        
    preds=pd.concat(preds,axis=0)
    preds['True positive']=preds['True positive'].fillna(False)
    preds['True negative']=preds['True negative'].fillna(False)
    preds['False positive']=preds['False positive'].fillna(False)
    preds['False negative']=preds['False negative'].fillna(False)
    
    qc=getModelQC(preds)
    print('R2 = '+str(qc['R2']))
    preds.to_csv(a.o+'/predictions_loo.csv')
    
    ########### Train the model on total set and export #############
    clf = svm.SVC(kernel='linear',class_weight=weights)
    clf.fit(data.drop(columns=['Group']), data['Group'])
    joblib.dump(clf, a.o+'/classifier.sav')
    
    ########### Train antimicrobial activity predictor independently #############
    data_balanced=pd.concat([data,data[data['Group']==1]])    #manually balancing the dataset by duplicating antimicrobial samples --> for the probability estimator ONLY
    probEst = svm.SVC(kernel='linear',class_weight={0:1,1:1})
    probEst.fit(data_balanced.drop(columns=['Group']), data_balanced['Group'])
    probEst= CalibratedClassifierCV(probEst, method='sigmoid', cv='prefit') # Platt scaling
    probEst.fit(data_balanced.drop(columns=['Group']), data_balanced['Group'])
    joblib.dump(probEst, a.o+'/probability_estimator.sav')
    

    ########### Export information about the predictor #############
    info_df=pd.DataFrame(index=['R2','parameters'],columns=['SVM'])
    info_df.loc['R2','SVM']=qc['R2']
    info_df.loc['Precision','SVM']=qc['Precision']
    info_df.loc['Recall','SVM']=qc['Recall']
    info_df.loc['Specificity','SVM']=qc['Specificity']
    info_df.loc['F-score','SVM']=qc['F1']
    info_df.loc['Balanced accuracy','SVM']=qc['Balanced R2']
    info_df.loc['Positive training set','SVM']=a.posPrefix
    info_df.loc['Negative training set','SVM']=a.negPrefix
    info_df.loc['Number of proteins in positive training set']=len(allProp[allProp['Group']==1])
    info_df.loc['Number of proteins in negative training set']=len(allProp[allProp['Group']==0])
    info_df.loc['Number of features used in training']=len(allProp.columns)-1
    if not a.excludeKmers:
        info_df.loc['K-mers used as features?']='Yes'
        info_df.loc['Number of k-mers used as features in training']=len(kept)
        info_df.loc['K-mers used as features in training']=';'.join(kept)
    else:
        info_df.loc['K-mers used as features?']='No'
    info_df.to_csv(a.o+'/classifier.info.csv')
	
    importances=clf.coef_[0]
    svm_importances = pd.Series(importances, index=data.drop(columns=['Group']).columns)
    svm_importances=pd.DataFrame(svm_importances)
    svm_importances=svm_importances.sort_values(by=svm_importances.columns[0],ascending=False)
    svm_importances.to_csv(a.o+'/classifier.feature_importances.csv')




