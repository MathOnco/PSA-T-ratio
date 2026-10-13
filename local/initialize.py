import pandas as pd
import numpy as np
from lmfit import Parameters

###############
# Load specs
###############
def getFittingSpecs(FittingSpecs):
    fitting0=FittingSpecs.DataFit.tolist()
    fitting0 = [x for x in fitting0 if str(x) != 'nan']
    fitting1=FittingSpecs.ModelFit.tolist()
    fitting1 = [x for x in fitting1 if str(x) != 'nan']
    fitting2=FittingSpecs.Weight.tolist()
    fitting2 = [x for x in fitting2 if str(x) != 'nan']
    fitting3=FittingSpecs.DetectionLimit.tolist()
    fitting3 = [x for x in fitting3 if str(x) != 'nan']
    fitting=[fitting0,fitting1,fitting2,fitting3]
    return [fitting]

def getPlotSpecs(modelSpecs):
    plotsIndex=modelSpecs.PlotsIndex.tolist()
    plotsIndex = [x for x in plotsIndex if str(x) != 'nan']
    plots=modelSpecs.Plots.tolist()
    plots = [x for x in plots if str(x) != 'nan']
    plotsData=modelSpecs.PlotsData.tolist()
    plotsData = [x for x in plotsData if str(x) != 'nan']
    plotsScale=modelSpecs.PlotsScale.tolist()
    plotsScale = [x for x in plotsScale if str(x) != 'nan']
    plotsTransform=modelSpecs.PlotsTransform.tolist()
    plotsTransform = [x for x in plotsTransform if str(x) != 'nan']
    return [plots,plotsData,plotsIndex,plotsScale,plotsTransform]

def getModelSpecs(modelSpecs):
    stateVarsNames=modelSpecs.StateVars.tolist()
    stateVarsNames = [x for x in stateVarsNames if str(x) != 'nan']
    stateVarsType=modelSpecs.StateVarsType.tolist()
    stateVarsType = [x for x in stateVarsType if str(x) != 'nan']
    pops0Input=modelSpecs.Pops0Input.tolist()
    pops0Input = [x for x in pops0Input if str(x) != 'nan']
    pops0Output=modelSpecs.Pops0Output.tolist()
    pops0Output = [x for x in pops0Output if str(x) != 'nan']
    pops0Type=modelSpecs.Pops0Type.tolist()
    pops0Type = [x for x in pops0Type if str(x) != 'nan']
    return [stateVarsNames,stateVarsType,pops0Input,pops0Output,pops0Type]

class Fitting2():
    def __init__(self,nIterations,n_sample,n_wins,fitting, **kwargs):
        self.nIterations=nIterations
        self.n_sample=n_sample
        self.n_wins=n_wins
        self.fitting=fitting

def defineInitValues(patientIdDic,dataDir,txCols,my_fit):
    tmp=[]
    for p in patientIdDic.keys():
        patientId = patientIdDic[p]
        patientDataDf = loadData(dataDir,patientId,txCols)
        T0=patientDataDf['Testosterone'].iloc[0]
        PSA0=patientDataDf['PSA'].iloc[0]

        tmp0=[]
        if('Testosterone' in my_fit.fitting[0]):
            tmp0=tmp0+[T0]
        if('PSA' in my_fit.fitting[0]):
            tmp0=tmp0+[PSA0]
        tmp.append(tmp0)

    # # list columns/names
    cols=[]
    if('Testosterone' in my_fit.fitting[0]):
        cols=cols + ['T0']
    if('PSA' in my_fit.fitting[0]):
        cols=cols + ['P0']

    lArray = pd.DataFrame(tmp,columns=cols)
    return lArray

def loadData(pathToData,patientId,txCols):
    patientDataDf = pd.read_excel(pathToData, sheet_name=str(patientId),header=1)
    dropCols=['Notes']
    for col in dropCols:
        patientDataDf.drop(col, axis=1, inplace=True)

    patientDataDfC=patientDataDf[patientDataDf['Cycle']>0]
    patientDataDfC = patientDataDfC.reset_index(drop=True)

    #get missing data, pretreatment values
    if(0 in patientDataDf['Cycle'].values):
        dataDf0=patientDataDf.copy()
        dataDf0=dataDf0[dataDf0['Cycle']==0]
        #Testo
        valsT0=[patientDataDfC['Testosterone'].iloc[0]]
        valsP0=[patientDataDfC['PSA'].iloc[0]]
        for tp in range(0,len(dataDf0['Days'])):
            valsT0.append(dataDf0['Testosterone'].iloc[tp])
            if(dataDf0['Days'].iloc[tp]>=-60):
                valsP0.append(dataDf0['PSA'].iloc[tp])
        tmpT=np.nanmean(valsT0)
        tmpP=np.nanmean(valsP0)
    else:
        tmpT=patientDataDfC['Testosterone'].iloc[0]
        tmpP=patientDataDfC['PSA'].iloc[0]
    T0=tmpT
    patientDataDfC.loc[0,'Testosterone']=T0
    PSA0=tmpP
    patientDataDfC.loc[0,'PSA']=PSA0

    #define drug index values
    strAll=[]
    for i in range(0,len(patientDataDfC[txCols[0]])):
        strNew=''
        for tx in txCols:
            strNew=strNew+str(patientDataDfC[tx].iloc[i])
        strAll.append(int(strNew,2))
    patientDataDfC['DrugIndex']=strAll
    
    return patientDataDfC

def loadDictionary(filename):
    fullDictionary= pd.read_excel(filename)
    dictionary = Parameters()
    for i in range(0,len(fullDictionary)):
        dictionary.add(fullDictionary.Name[i], value=fullDictionary.Value[i], min=fullDictionary.Min[i], max=fullDictionary.Max[i], vary = True)
    print("done: initialize parameters")
    return dictionary

def loadFullDictionary(modelName): 
    dicFile = './Parameters/Pars_'+modelName+'.xlsx'
    fullDictionary= pd.read_excel(dicFile,'Pars')
    parDic = Parameters()
    for i in range(0,len(fullDictionary)):
        if(fullDictionary.Type[i]=='C'):
            parDic.add(fullDictionary.Name[i], value=fullDictionary.Value[i], min=fullDictionary.Min[i], max=fullDictionary.Max[i], vary = False)
        else:
            parDic.add(fullDictionary.Name[i], value=fullDictionary.Value[i], min=fullDictionary.Min[i], max=fullDictionary.Max[i], vary = True)
    return [fullDictionary,parDic]


