import numpy as np
import pandas as pd
import scipy
import sys 
from lmfit import minimize, Parameters

sys.path.append("./local")
import myUtils as utils
import odeFuns as odeFuns
    
def runSim(tmpModel,patientId,patInd,lArray,parDicAdd,tx,my_fit,boolTinput,testo,solver_kws):
    print(patientId)
    params = Parameters() #initialize parameters
    params=utils.setConParamsFromDictionary(tmpModel.varsC,params,tmpModel.dictionary)
    for par in lArray.columns:
        params.add(par, value=lArray[par][patInd],vary=False)
    for par in parDicAdd.keys():
        params.add(par, value=parDicAdd[par],vary=False)

    #get and set population splits based on burden
    [consSet,consSplit,outsSplit] = utils.getSplits(tmpModel)
    params = utils.setSplits(consSplit,lArray,patInd,params,outsSplit)
    tmpModel.SetParams(**params.valuesdict())

    if(boolTinput):
        tmpModel.SimulateWithT(testo,treatmentScheduleList=tx,**solver_kws)
    else:
        tmpModel.Simulate(treatmentScheduleList=tx,**solver_kws)
    topDf = tmpModel.resultsDf
    return topDf

################
# Data analysis
################
def find0(feature,dataDf,colsDf,valsDf):
    nameDf=feature+'0'
    colsDf.append(nameDf)
    f0=dataDf[feature].loc[0]
    valsDf.append(f0)
    return [f0,colsDf,valsDf]

def findMean(feature,dataDf,colsDf,valsDf):
    nameDf=feature+'Mean'
    colsDf.append(nameDf)
    valsDf.append(np.nanmean(dataDf[feature]))
    return [colsDf,valsDf]

def findMax(feature,dataDf,colsDf,valsDf):
    nameDf=feature+'Max'
    colsDf.append(nameDf)
    valsDf.append(np.nanmax(dataDf[feature]))
    return [colsDf,valsDf]

def findMin(feature,dataDf,colsDf,valsDf):
    nameDf=feature+'Min'
    colsDf.append(nameDf)
    valsDf.append(np.nanmin(dataDf[feature]))
    return [colsDf,valsDf]

def findSwitchIndex(feature,dataDf):
    if(0 in dataDf['DrugIndex'].values and len(dataDf['Days'])>1):
        indStmp=np.where(dataDf['DrugIndex']==0)[0]
        indSwitch=indStmp[0]
        if(np.isnan(dataDf[feature].iloc[indSwitch])):
            dataDf['nearSwitch']=abs(dataDf['Days']-dataDf['Days'].iloc[indSwitch])
            sorted_df = dataDf.sort_values(by='nearSwitch', ascending=True)
            vi=0
            while(np.isnan(sorted_df[feature].iloc[vi])):
                vi=vi+1            
            indSwitch=sorted_df.index[vi]
    return indSwitch

def findSwitchVal(feature,dataDf,ind,colsDf,valsDf):
    name=feature+'sw'
    fS=dataDf['PSA'].iloc[ind]
    colsDf.append(name)
    valsDf.append(fS)
    return [fS,colsDf,valsDf]

def findFinalVal(feature,dataDf,colsDf,valsDf):
    name=feature+'Final'
    fF=dataDf[feature].iloc[-1]
    if(np.isnan(fF)):
        fF=dataDf[feature].iloc[-2]
    colsDf.append(name)
    valsDf.append(fF)
    return [fF,colsDf,valsDf]   

def calculateMetrics(pat,pi,dataDf,dataVals,modVals,fitVals):
    colsDf=[]
    valsDf=[]
    [PSA0,colsDf,valsDf]=find0('PSA',dataDf,colsDf,valsDf)
    [T0,colsDf,valsDf]=find0('Testosterone',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMean('PSA',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMean('Testosterone',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMax('PSA',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMax('Testosterone',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMin('PSA',dataDf,colsDf,valsDf)
    [colsDf,valsDf]=findMin('Testosterone',dataDf,colsDf,valsDf)

    indSwitch=findSwitchIndex('PSA',dataDf)
    [PSAS,colsDf,valsDf]=findSwitchVal('PSA',dataDf,indSwitch,colsDf,valsDf)
    [TS,colsDf,valsDf]=findSwitchVal('Testosterone',dataDf,indSwitch,colsDf,valsDf)
    [PSAF,colsDf,valsDf]=findFinalVal('PSA',dataDf,colsDf,valsDf)
    [TF,colsDf,valsDf]=findFinalVal('Testosterone',dataDf,colsDf,valsDf)

    #cycle data
    DaysOn=dataDf['Days'].iloc[indSwitch]/30
    DaysOff=(max(dataDf['Days']) - dataDf['Days'].iloc[indSwitch-1])/30 
    DaysTotal=max(dataDf['Days'])/30
    colsDf.append('DaysOn')
    valsDf.append(DaysOn)
    colsDf.append('OffTime')
    valsDf.append(DaysOff)#actually in months
    colsDf.append('On/Off')
    valsDf.append(DaysOn/DaysOff)
    colsDf.append('CycleTime')
    valsDf.append(DaysTotal)# in months
    
    #PSA metrics
    #drug on
    colsDf.append('diffPOnRate')
    valsDf.append((PSAS-PSA0)/(DaysOn))
    colsDf.append('diffPOn')
    valsDf.append(PSAS-PSA0)
    colsDf.append('diffPOnRatePer')
    valsDf.append((PSAS-PSA0)/(PSA0*DaysOn))
    colsDf.append('diffPOnPer')
    valsDf.append((PSAS-PSA0)/PSA0)
    colsDf.append('diffPOnPerPos')
    valsDf.append(100*np.abs(PSAS-PSA0)/PSA0)
    #drug off
    colsDf.append('diffPOffRate')
    valsDf.append((PSAF-PSAS)/DaysOff)
    colsDf.append('diffPOff')
    valsDf.append(PSAF-PSAS)
    colsDf.append('diffPOffRatePer')
    valsDf.append((PSAF-PSAS)/(PSAS*DaysOff))
    colsDf.append('diffPOffPer')
    valsDf.append((PSAF-PSAS)/PSA0)
    #total cycle
    colsDf.append('diffPTotRate')
    valsDf.append((PSAF-PSA0)/DaysTotal)
    colsDf.append('diffPTot')
    valsDf.append(PSAF-PSA0)
    colsDf.append('diffPTotRatePer')
    valsDf.append((PSAF-PSA0)/(PSA0*DaysTotal))
    colsDf.append('diffPTotPer')
    valsDf.append((PSAF-PSA0)/PSA0)

    #testosterone metrics
    #drug on
    colsDf.append('diffTOnRate')
    valsDf.append((TS-T0)/(DaysOn))
    colsDf.append('diffTOn')
    valsDf.append((TS-T0))
    colsDf.append('diffTOnRatePer')
    valsDf.append((TS-T0)/(DaysOn*T0))
    colsDf.append('diffTOnPer')
    valsDf.append((TS-T0)/T0)
    #drug off
    colsDf.append('diffTOffRate')
    valsDf.append((TF-TS)/(DaysOff))
    colsDf.append('diffTOff')
    valsDf.append((TF-TS))
    colsDf.append('diffTOffRatePer')
    valsDf.append((TF-TS)/(TS*DaysOff))
    colsDf.append('diffTOffPer')
    valsDf.append((TF-TS)/TS)
    #total
    colsDf.append('diffTTotRate')
    valsDf.append((TF-T0)/(DaysTotal))
    colsDf.append('diffTTot')
    valsDf.append((TF-T0))
    colsDf.append('diffTTotRatePer')
    valsDf.append((TF-T0)/(T0*DaysTotal))
    colsDf.append('diffTTotRatePerPos')
    valsDf.append(100*np.abs(TF-T0)/(T0*DaysTotal))
    colsDf.append('diffTTotPer')
    valsDf.append((TF-T0)/T0)
    colsDf.append('diffTTotPerPos')
    valsDf.append(100*np.abs(TF-T0)/T0)

    #interpolated
    teval=range(0,max(dataDf['Days']))
    nanT=np.isnan(dataDf['Testosterone'])
    if(len(dataDf.Days[~nanT])>1):
        f1T = scipy.interpolate.interp1d(dataDf.Days[~nanT],dataDf['Testosterone'][~nanT],fill_value="extrapolate",kind="linear")#"extrapolate"
        tmpTest=f1T(teval)
        daysThigh=0
        for day in teval:
            if(tmpTest[day]>=100):daysThigh+=1
        colsDf.append('PerTimeTHighTotal')
        valsDf.append(100*daysThigh/max(dataDf['Days']))
        

    #P/T
    fact=.01
    colsDf.append('PSA0/T0')
    valsDf.append(fact*PSA0/T0)
    colsDf.append('PSAsw/Tsw')
    valsDf.append(fact*PSAS/TS)
    colsDf.append('PSAF/TF')
    valsDf.append(fact*PSAF/TF)
    colsDf.append('avePSA/T')
    valsDf.append(fact*np.nanmean(dataDf['PSA']/dataDf['Testosterone']))## back to interpolated values
    colsDf.append('maxPSA/T')
    valsDf.append(fact*np.nanmax(dataDf['PSA']/dataDf['Testosterone']))## back to interpolated values
    #drug on
    colsDf.append('diffOnRate')
    valsDf.append(fact*(PSAS/TS-PSA0/T0)/(DaysOn))
    colsDf.append('diffOn')
    valsDf.append(fact*(PSAS/TS-PSA0/T0))
    colsDf.append('diffOnRatePer')
    valsDf.append(fact*(PSAS/TS-PSA0/T0)/((PSA0/T0)*DaysOn))
    colsDf.append('diffOnPer')
    valsDf.append(fact*(PSAS/TS-PSA0/T0)/((PSA0/T0)))
    colsDf.append('diffOnPer2')
    valsDf.append(100*fact*(PSAS/TS-PSA0/T0)/((PSA0/T0)))
    #drug off
    colsDf.append('diffOffRate')
    valsDf.append(fact*(PSAF/TF-PSAS/TS)/(DaysOff))
    colsDf.append('diffOff')
    valsDf.append(fact*(PSAF/TF-PSAS/TS))
    colsDf.append('diffOffRatePer')
    valsDf.append(fact*(PSAF/TF-PSAS/TS)/((PSAS/TS)*DaysOff))
    colsDf.append('diffOffPer')
    valsDf.append(fact*(PSAF/TF-PSAS/TS)/((PSAS/TS)))
    #Total
    colsDf.append('diffTotRate')
    valsDf.append(fact*(PSAF/TF-PSA0/T0)/(DaysTotal))
    colsDf.append('diffTot')
    valsDf.append(fact*(PSAF/TF-PSA0/T0))
    colsDf.append('diffTotRatePer')
    valsDf.append(fact*(PSAF/TF-PSA0/T0)/((PSA0/T0)*DaysTotal))
    colsDf.append('diffTotPer')
    valsDf.append(fact*(PSAF/TF-PSA0/T0)/((PSA0/T0)))

    colsDf.append('Gleason')
    valsDf.append(dataVals['GS'][pi])
    colsDf.append('BMs')
    valsDf.append(dataVals['BMs'][pi])
    colsDf.append('LNs')
    valsDf.append(dataVals['LNs'][pi])
    colsDf.append('Total mets')
    valsDf.append(dataVals['TMs'][pi])
    colsDf.append('risk')
    valsDf.append(dataVals['risk'][pi])
    colsDf.append('de novo')
    valsDf.append(dataVals['dNM1'][pi])
    colsDf.append('Age')
    valsDf.append(dataVals['Age'][pi])
    colsDf.append('BMI')
    valsDf.append(dataVals['BMI'][pi])

    colsDf.append('maxTfit')
    valsDf.append(fitVals['maxTfit'][pi])
    colsDf.append('gRTfit')
    valsDf.append(fitVals['gRTfit'][pi])
    colsDf.append('tdelayfit')
    valsDf.append(fitVals['tdelayfit'][pi])
    colsDf.append('maxPfit')
    valsDf.append(fitVals['maxPfit'][pi])
    colsDf.append('dRPfit')
    valsDf.append(fitVals['dRPfit'][pi])
    colsDf.append('minTfit')
    valsDf.append(fitVals['minTfit'][pi]) 
    colsDf.append('dRTfit')
    valsDf.append(fitVals['dRTfit'][pi])
    colsDf.append('minPfit')
    valsDf.append(fitVals['minPfit'][pi])
    colsDf.append('gRPfit')
    valsDf.append(fitVals['gRPfit'][pi])

    colsDf.append('rS')
    valsDf.append(modVals['rS'][pi])
    colsDf.append('X0')
    valsDf.append(modVals['X0'][pi])

    colsDf.append('id')
    valsDf.append(pat)
    colsDf.append('progPSA')
    valsDf.append(dataVals['progPSA'][pi])
    colsDf.append('progRP')
    valsDf.append(dataVals['progRP'][pi])
    colsDf.append('progOS')
    valsDf.append(dataVals['progOS'][pi])
    colsDf.append('dAbi')
    valsDf.append(dataVals['Abi'][pi])
    colsDf.append('TTP')
    valsDf.append(dataVals['TTP'][pi])

    pData=pd.DataFrame()
    pData['names']=colsDf
    pData['values']=valsDf
    pData['id']=pat
    pData['ELN']=dataVals['ELNstatus'][pi]
    pData['OS']=dataVals['Dstatus'][pi]
    pData['PN']=dataVals['PNstatus'][pi]
    pData['RP']=dataVals['RPstatus'][pi]
    pData['dAbi']=dataVals['Abi'][pi]
    return pData


def calculateMetrics2(pat,pi,dataDf,dataVals,modVals,fitVals):
    colsDf=[]
    valsDf=[]
    PSA0=dataDf['PSA'].loc[0]
    colsDf.append('PSA0')
    valsDf.append(PSA0)
    colsDf.append('aveP')
    valsDf.append(np.nanmean(dataDf['PSA']))

    tmpD=dataDf['Days'][len(dataDf['Days'])-1]
    indF=tmpD-1
    teval=range(0,tmpD)
    nanP=np.isnan(dataDf['PSA'])
    if(len(dataDf.Days[~nanP])>1):
        f1P = scipy.interpolate.interp1d(dataDf.Days[~nanP],dataDf['PSA'][~nanP],fill_value="extrapolate",kind="linear")
        tmpPSA=f1P(teval)
    nanT=np.isnan(dataDf['Testosterone'])
    if(len(dataDf.Days[~nanT])>1):
        f1T = scipy.interpolate.interp1d(dataDf.Days[~nanT],dataDf['Testosterone'][~nanT],fill_value="extrapolate",kind="linear")
        tmpTest=f1T(teval)

    #calculated time point indexes
    indFinal=len(dataDf['DrugIndex'])-1
    if(0 in dataDf['DrugIndex'].values and len(dataDf['Days'])>1):
        indStmp=np.where(dataDf['DrugIndex']==0)[0]
        if(len(indStmp)>1):
            indSwitch=indStmp[0]
            indS=dataDf['Days'].iloc[indSwitch]
        elif(len(indStmp)==1):
            if isinstance(indStmp, np.ndarray):
                indSwitch=indStmp[0]
            else:
                indSwitch=indStmp
            indS=dataDf['Days'].iloc[indSwitch]
        else:
            indSwitch=0
            indS=0
    else:
        indSwitch=0
        indS=0

    pData=pd.DataFrame()
    pData['names']=colsDf
    pData['values']=valsDf
    pData['id']=pat
    pData['progressed1']=dataVals['prog1'][pi]
    pData['progressed2']=dataVals['prog2'][pi]
    pData['progressed3']=dataVals['prog3'][pi]
    pData['progressed4']=dataVals['prog4'][pi]
    pData['dAbi']=dataVals['drug'][pi]
    return pData
  
print("done import aFuns")
