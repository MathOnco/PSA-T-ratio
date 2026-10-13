#myPlots.py
import os
import matplotlib as mpl
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
plt.rcParams.update({
    'axes.titlesize': 26,     # Title
    'axes.labelsize': 24,     # X and Y labels
    'xtick.labelsize': 22,    # X tick labels
    'ytick.labelsize': 22,    # Y tick labels
    'legend.fontsize': 14,     # Legend
    'xtick.major.pad':1,
    'ytick.major.pad':1
})
mpl.rcParams['axes.formatter.limits'] = (-99, 99)
from lmfit import minimize, Parameters
import scipy
import seaborn as sns   #plot options - ease and beauty 
import pandas as pd  #data frames
import numpy as np   #arrays/matrixes - matlab
import math as math
import myUtils as utils
import initialize as initialize


def plot_sims(parsDf,lArray,index_pat, myModel, modelSpecs, my_fit, dataDf,boolTinput,fitting,tx,solver_kws={}, ax=None, **kwargs):
    print("plotting....")
    dataFit=fitting[0]
    modelFit=fitting[1]
    fitWeight=fitting[2]

    [plots,plotsData,plotsIndex,plotsScale,plotsTransform] = initialize.getPlotSpecs(modelSpecs)

       
    #setup plots
    popCount = 1 if (myModel.stateVarsTypeList.count("Pop")>=1 or myModel.stateVarsTypeList.count("Pop1")>=1) else 0
    otherCount = len(np.unique(plotsIndex))-popCount
    numPlots = 1 + popCount + otherCount
    totalHeight = 1 + 3*popCount + 3*otherCount
    height_ratios = [.5]
    for ps in range(0,popCount):
        height_ratios.append(1.5)
    for ps in range(0,otherCount):
        height_ratios.append(1.5)
    fig, ax = plt.subplots(numPlots, 1, figsize=kwargs.get("figsize",(10,totalHeight)), gridspec_kw={'height_ratios': height_ratios})#, sharex=True
    fig.subplots_adjust(hspace=0)
    
    # Simulate top 10
    topDf=[]
    topDf=pd.DataFrame(topDf)
    minT=[]
    maxT=[]
    params = Parameters() 
    [consSet,consSplit,outsSplit] = utils.getSplits(myModel)

    for i in range(0,parsDf.shape[0]):
        params=utils.setParams(params,myModel.varsC,myModel.varsV,parsDf,myModel.dictionary,i) ##clean up        
        params = utils.setConstants(consSet,lArray,params,index_pat,myModel.dictionary)
        params = utils.setSplits(consSplit,lArray,index_pat,params,outsSplit)

        myModel.SetParams(**params.valuesdict())
        if(boolTinput):
            nanT=np.isnan(dataDf['Testosterone'])
            dataDf['Tlog']=np.log10(dataDf['Testosterone'])
            f1T = scipy.interpolate.interp1d(dataDf.Days[~nanT],dataDf['Tlog'][~nanT],fill_value="extrapolate",kind="linear")#"extrapolate"
            teval=range(0,dataDf['Days'][len(dataDf['Days'])-1])
            testo=f1T(teval)
            testo=10**testo
            i=0
            tx=[]
            for d in teval:
                tx0=[]
                tx0.append(d)
                tx0.append(d+1)
                for col in myModel.txCols:
                    tx0.append(dataDf[col].loc[i])
                tx.append(tx0)
                if(d>=dataDf['Days'].loc[i+1]):
                    i=i+1 

            myModel.SimulateWithT(testo,treatmentScheduleList=tx,max_step=1,**solver_kws)
        else:
            myModel.Simulate(treatmentScheduleList=tx,max_step=1,**solver_kws)
        myModel.resultsDf['id']=i
        topDf = pd.concat([topDf, myModel.resultsDf],ignore_index=True)
        l1=myModel.resultsDf[modelFit[0]].values
        l1=l1[l1 < 1E308]
        l1=l1[l1 > -1E308]
        minT.append(np.min(l1))
        maxT.append(np.max(l1))
    xmin=0
    xmax=max(myModel.resultsDf['Days'])
    xmarg=1 

    #plotting
    axi=0
    axx=ax[axi]
    plot_dose_over_time(axx,myModel,tx,dataDf,xmin,xmax,xmarg)
    axx.tick_params(labelbottom=False)
    
    for j in range(0,numPlots-1):
        axi=axi+1
        axx=ax[axi]
        p = [i for i, x in enumerate(plotsIndex) if x==j]
        for pi in p:
            if(plots[pi]=='TumorSize'):
                plot_tumor_over_time(axx,xmin,xmax,topDf,kwargs,dataDf,xmarg,myModel,modelFit,dataFit,plots,plotsData,plotsScale,plotsIndex,pi)
            else:
                metric=plots[pi]
                plot_metric_over_time(axx,xmin,xmax,topDf,kwargs,dataDf,xmarg,metric,myModel,modelFit,dataFit,plots,plotsData,plotsScale,plotsTransform)
        if(j<numPlots-2):axx.tick_params(labelbottom=False)          
        
####################
# Plot types 
####################

def plot_dose_over_time(axx,myModel,tx,dataDf,xmin,xmax,xmarg):
    cols=['ti','tf']+myModel.txCols
    tx=pd.DataFrame(tx,columns=cols)
    for i,drug in enumerate(myModel.txCols):
        drugBarPosition = 1-.5*i
        plot_drug(tx,drug,drugBarPosition,axx,color=myModel.colorPalette[drug])
    axx.set_xlim((xmin-xmarg, xmax+xmarg))
    axx.tick_params(labelleft=False) 

def plot_drug(tx,drug,drugBarPosition,ax,color):  
    drugConcentrationVec=tx[drug]
    drugConcentrationVec = np.array([x/(np.max(drugConcentrationVec)+1e-12) for x in drugConcentrationVec])
    drugConcentrationVec = 0.5*drugConcentrationVec + drugBarPosition
    axName=ax
    axName.fill_between(tx['ti'], drugConcentrationVec, drugBarPosition, 
                        step="post", color=color, alpha=0.6, label=drug)

def plot_tumor_over_time(axx,xmin,xmax,topDf,kwargs,dataDf,xmarg,myModel,modelFit,dataFit,plots,plotsData,plotsScale,plotsIndex,pi):
    popName=plots[pi]
    plotMods(topDf,['Days'],[popName],5,axx,myModel.colorPalette,0.5,**kwargs)
    tumorPopsList = []
    tmp=popName.replace('TumorSize','')
    popName='Pop'+tmp
    sizeName='TumorSize'+tmp
    for p,pops in enumerate(myModel.stateVarsTypeList):
        if(pops == 'Pop' or pops==popName):
            tumorPopsList.append(myModel.stateVars[p])
    for pop in tumorPopsList:
        plotMods(topDf,['Days'],[pop],3,axx,myModel.colorPalette,0.5,**kwargs)

    ind2=plots.index(sizeName)
    if(plotsData[ind2]!='Non'):
        sns.scatterplot(data=dataDf,x='Days', y=plotsData[ind2], color=myModel.colorPalette[sizeName],edgecolor="black",ax=axx)
    plt.ylabel("Tumor size") 
    axx.set_xlim((xmin-xmarg, xmax+xmarg))
    axx.get_legend().remove()
    if(plotsScale[ind2]=='log'):
        axx.set_yscale('log',base=10)
        formatter = ticker.ScalarFormatter()
        formatter.set_scientific(False)
        axx.yaxis.set_major_formatter(formatter)
        axx.yaxis.set_minor_formatter(formatter)
        axx.yaxis.set_minor_locator(ticker.MaxNLocator(nbins=3)) 

def plotMods(df,tCol,modFits,lw0,ax0,palette,alpha):
    modelFitsDf = pd.melt(df, id_vars=tCol, value_vars=modFits)
    sns.lineplot(x=tCol[0],y="value", hue="variable", style="variable",lw=lw0, palette=palette, alpha=alpha,legend='full', data=modelFitsDf, ax=ax0)

def plot_metric_over_time(axx,xmin,xmax,topDf,kwargs,dataDf,xmarg,metric,myModel,modelFit,dataFit,plots,plotsData,plotsScale,plotsTransform):
    ind2=plots.index(metric)
    if(plotsTransform[ind2]=="norm" and (metric in modelFit)):
        ind = modelFit.index(metric)
        normName=metric+"norm"  
        data1=dataDf[dataFit[ind]].interpolate(method ='linear', limit_direction ='backward')  
        topDf[normName] = topDf[metric]*data1.iloc[0]/topDf[metric][0]
        plotMods(topDf,['Days'],[normName],5,axx,myModel.colorPalette,0.5,**kwargs)
    elif(plotsTransform[ind2]=="splitnorm" and (metric in modelFit)):
        ind = modelFit.index(metric)
        normName="splitnorm"  
        data1=dataDf[dataFit[ind]].interpolate(method ='linear', limit_direction ='backward')  
        topDf[normName] = topDf['split']*data1.iloc[0]/topDf['split'][0]
        plotMods(topDf,['Days'],[normName],5,axx,myModel.colorPalette,0.5,**kwargs)
    else:
        plotMods(topDf,['Days'],[metric],5,axx,myModel.colorPalette,0.5,**kwargs)
    if(plotsData[ind2]!='None'):
        sns.scatterplot(data=dataDf,x='Days', y=plotsData[ind2], color=myModel.colorPalette[plotsData[ind2]],edgecolor="black",ax=axx)
    plt.ylabel(metric) 
    axx.set_xlim((xmin-xmarg, xmax+xmarg))
    if(plotsScale[ind2]=='log'):
        axx.set_yscale('log',base=10)
        formatter = ticker.ScalarFormatter()
        formatter.set_scientific(False)
        axx.yaxis.set_major_formatter(formatter)
        axx.yaxis.set_minor_formatter(formatter)
        axx.yaxis.set_minor_locator(ticker.MaxNLocator(nbins=3)) 

def plot_pars(patsdf,lPars,dictionary):
    ncolMax=6
    colSize=2
    if(len(lPars)<ncolMax):
        ncol=len(lPars)
    else:
        ncol=ncolMax
    nrow=math.floor((len(lPars)-1)/ncol)+1

    x1=colSize*ncolMax*2
    x2=colSize*nrow
    fig, axis = plt.subplots(nrow, ncol,figsize=(x1,x2))

    ax1=0
    ax2=0
    for i in range(0,len(lPars)): 
        ff=dictionary.get(lPars[i])
        if(nrow==1):
            axis[ax2].hist(patsdf[lPars[i]], bins=10) 
            axis[ax2].set_title(lPars[i])
        else:  
            axis[ax1, ax2].hist(patsdf[lPars[i]], bins=10)
            axis[ax1, ax2].set_title(lPars[i])
        ax2=ax2+1
        if(ax2>=ncolMax):
            ax2=0
            ax1=ax1+1
        
def plot_pars_over_time(pars,path,dictionary,nIterations):
    meanTime=[]
    stdTime=[]
    xTime=[]
    for i in range(1,nIterations+1):
        dd=path + str(i) + '.xlsx'
        parsI = pd.read_excel(dd,index_col=[0])
        tmp = parsI.mean(axis=0)
        tmp2 = parsI.std(axis=0)
        meanTime.append(tmp.values.tolist())
        stdTime.append(tmp2.values.tolist())
        xTime.append(i)
    meanTime=pd.DataFrame(meanTime)
    stdTime=pd.DataFrame(stdTime)

    ncolMax=6
    colSize=2
    if(len(pars)<ncolMax):
        ncol=len(pars)
    else:
        ncol=ncolMax
    nrow=math.floor((len(pars)-1)/ncol)+1

    x1=colSize*ncolMax*2
    x2=colSize*nrow
    fig, axis = plt.subplots(nrow, ncol,figsize=(x1,x2))

    ax1=0
    ax2=0
    for i, par in enumerate(pars):
        if(nrow==1):
            axis[ax2].fill_between(xTime,meanTime[i]-stdTime[i], meanTime[i]+stdTime[i],alpha=0.5)
            axis[ax2].plot(xTime,meanTime[i])
            axis[ax2].set_title(pars[i])
        else: 
            axis[ax1, ax2].fill_between(xTime,meanTime[i]-stdTime[i], meanTime[i]+stdTime[i],alpha=0.5) 
            axis[ax1, ax2].plot(xTime,meanTime[i])        
            axis[ax1, ax2].set_title(pars[i])
        ax2=ax2+1
        if(ax2>=ncolMax):
            ax2=0
            ax1=ax1+1


def plot_pars_dists(pars,path,dictionary,nIterations,lesionId,path_to_file,n_wins): 
    dd=path_to_file
    parsI = pd.read_excel(dd,index_col=[0])

    ncolMax=6
    colSize=2
    if(len(pars)<ncolMax):
        ncol=len(pars)
    else:
        ncol=ncolMax
    nrow=math.floor((len(pars)-1)/ncol)+1

    x1=colSize*ncolMax*2
    x2=colSize*nrow
    fig, axis = plt.subplots(nrow, ncol,figsize=(x1,x2))

    ax1=0
    ax2=0
    for par in pars:
        if(len(pars)==1):
            axis.set_title(par)
            axis.hist(parsI[par][:n_wins].values,range=[dictionary.get(par).min,dictionary.get(par).max])
        elif(nrow==1):
            axis[ax2].set_title(par)
            axis[ax2].hist(parsI[par][:n_wins].values,range=[dictionary.get(par).min,dictionary.get(par).max])
        else: 
            axis[ax1, ax2].set_title(par)
            axis[ax1, ax2].hist(parsI[par][:n_wins].values,range=[dictionary.get(par).min,dictionary.get(par).max])       
        ax2=ax2+1
        if(ax2>=ncolMax):
            ax2=0
            ax1=ax1+1

    pf_name="pars"+"%s.pdf"%lesionId
    path_to_fig=os.path.join(path,pf_name)    
    plt.savefig(path_to_fig)
    plt.close()