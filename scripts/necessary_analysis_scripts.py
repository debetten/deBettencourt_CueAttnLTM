import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import os, sys, pickle
from scipy.stats import shapiro, normaltest, anderson, ttest_1samp, mannwhitneyu, wilcoxon

def prettify_plot(ax,
					 xlim=None,xt=None,xtl=None,xl=None,xaxoffset=None,
					 ylim=None,yt=None,ytl=None,yl=None,ylrot=None,yaxoffset=None,
					 t=None,legend=None,legendloc=None):
    '''
    This is a plotting script that makes the default matplotlib plots a little prettier
    '''

    if os.path.isfile("/Library/Fonts/HelveticaNeue-Light.ttf"): 
        prop_light = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue-Light.ttf')    
    else: 
        prop_light = fm.FontProperties()

    if os.path.isfile("/Library/Fonts/HelveticaNeue.ttf"): 
        prop_reg = fm.FontProperties(fname='/Library/Fonts/HelveticaNeue.ttf')    
    else: 
        prop_reg = fm.FontProperties()

    ax.spines['bottom'].set_linewidth(1)
    ax.spines['bottom'].set_color("gray")
    if xaxoffset is not None: ax.spines['bottom'].set_position(('outward', 10))
    if yaxoffset is not None: ax.spines['left'].set_position(('outward', 10))

    ax.spines['left'].set_linewidth(1)
    ax.spines['left'].set_color("gray")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.yaxis.set_ticks_position("left")   
    ax.tick_params(axis='y',direction='out',length=5,width=1,color='gray')
    ax.xaxis.set_ticks_position("bottom") 
    ax.tick_params(axis='x',direction='out',length=5,width=1,color='gray')
    
    if yt is not None: ax.set_yticks(yt)
    if ytl is not None: ax.set_yticklabels((ytl),fontsize=32,fontproperties=prop_light) 
    if yl is not None: h = ax.set_ylabel(yl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if ylim is not None: ax.set_ylim(ylim) 
        
    if xt is not None: ax.set_xticks(xt)
    if xtl is not None: ax.set_xticklabels((xtl),fontsize=32,fontproperties=prop_light)
    if xl is not None: ax.set_xlabel(xl,fontsize=36,fontproperties=prop_reg,labelpad=12)
    if xlim is not None: ax.set_xlim(xlim)
    

    if t is not None: ax.set_title(t,y=1.08,fontsize=36,fontproperties=prop_reg)
    if legend is not None: 
        if legendloc is None: L = ax.legend(loc='center left', bbox_to_anchor=(0,.85))
        else: L = ax.legend(loc='center right', bbox_to_anchor=legendloc)
        plt.setp(L.texts,fontsize='large',fontproperties=prop)
    ax.tick_params(axis='both',pad=10)
    plt.locator_params(nbins=8)
    plt.tight_layout()

def plot_adjust_spines(ax, spines):
    '''
    This is another plotting script, that moves spines outward if you don't want the x and y axes to touch 
    '''
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points

@np.vectorize
def calculate_aprime(h,fa):
    '''
    Calculates sensitivity, according to the aprime formula
    '''
    if np.greater_equal(h,fa): a = .5 + (((h-fa) * (1+h-fa)) / (4 * h * (1-fa)))
    else: a = .5 - (((fa-h) * (1+fa-h)) / (4 * fa * (1-h)))
    return a

def resampling_statistics(data,chance,nsubj=None,nsamples=100000,skipfig=True):
    '''
    Nonparametric resampling statistics
    '''
    if nsubj is None: nsubj = np.size(data)
        
    #resample subjects with replacement 100,000 times
    subj_resampled = np.random.randint(0,nsubj,(nsamples,nsubj)) 

    #the empy matrix of resampled data for each of these resampled iterations
    data_resampled = np.empty(nsamples) 

    #recalculate mean given the resampled subjects
    for i in range(0,nsamples): data_resampled[i] = np.mean(data[subj_resampled[i,:]])

    #calculate p value 
    p = np.sum(np.less(data_resampled,chance))/float(nsamples) #count number of resample iterations below chance
    if np.equal(p,0): p = 1./float(nsamples)

    if skipfig is False:
        plt.figure(figsize=(4,3))
        ax = plt.subplot(111)
        plt.hist(data_resampled,normed=0,facecolor='gray',edgecolor='gray')
        plt.axvline(np.mean(data_resampled),color='b',lw=2,label='resampled mean')
        plt.axvline(np.mean(data),color='m',lw=1,label='original mean')
        plt.axvline(chance,color='c',lw=2,label='chance')
        make_plot_pretty(ax,ylrot=90,yl='Count (#)',legend='1') 
        plt.show()
    
    return p, data_resampled

def resampling_statistics_betweengroups(data1,data2,nsamples=100000,skipfig=True):
    '''
    Nonparametric resampling statistics - doesn't actually work yet
    '''
    nsubj1 = np.size(data1)
    nsubj2 = np.size(data2)
        
    #resample subjects with replacement 100,000 times within group
    subj1_resampled = np.random.randint(0,nsubj1,(nsamples,nsubj1)) 
    subj2_resampled = np.random.randint(0,nsubj2,(nsamples,nsubj2)) 

    #the empy matrix of resampled data for each of these resampled iterations
    data1_resampled = np.empty(nsamples) 
    data2_resampled = np.empty(nsamples) 

    #recalculate mean given the resampled subjects
    for i in range(0,nsamples): 
        data1_resampled[i] = np.mean(data1[subj1_resampled[i,:]])
        data2_resampled[i] = np.mean(data2[subj2_resampled[i,:]])

    #calculate p value 
    p = np.sum(np.less(data1_resampled,data2_resampled))/float(nsamples) #count number of resample iterations below chance
    if np.equal(p,0): p = 1./float(nsamples)

    if skipfig is False:
        plt.figure(figsize=(4,3))
        ax = plt.subplot(111)
        plt.hist(data_resampled,normed=0,facecolor='gray',edgecolor='gray')
        plt.axvline(np.mean(data_resampled),color='b',lw=2,label='resampled mean')
        plt.axvline(np.mean(data),color='m',lw=1,label='original mean')
        plt.axvline(chance,color='c',lw=2,label='chance')
        make_plot_pretty(ax,ylrot=90,yl='Count (#)',legend='1') 
        plt.show()
    
    return p, data1_resampled, data2_resampled

def load_data(project_name='sustAttnWM01',base_project_dir = '/Users/megan/Documents/projects/', pickle_files =True):
    '''
    Load experimental data
    '''
    
    #project directory
    project_dir = base_project_dir + project_name + '/' 

    #subject directory
    subjects_dir = project_dir + 'subjects/'
    subj_name = [f for f in os.listdir(subjects_dir) if f.endswith(project_name)]
    
    nsubj = np.size(subj_name)

    if pickle_files:
        display_dir = project_dir + 'display/'
        sys.path.insert(0, display_dir)
        import settings_sustAttnWM
    
        #fill in subj_dat data structure with the data from the saved pickle file for each subject
        subj_dat = {}
        for isubj in range(nsubj):
            subj_datadir = project_dir + 'subjects/' + subj_name[isubj] + '/data/beh/'
            subj_dat[isubj] = pickle.load( open(subj_datadir + subj_name[isubj] + '_expdat.p','rb'))

    else:
        pass 
        #need to write this function to load CSV files

    return subj_dat

def test_normality(data):
    #conducts 3 different tests of normality

    normal_test_results = [1,1,1]
    alpha = 0.05

    # Shapiro-Wilk Test
    stat, p = shapiro(data)
    if p > alpha:
        print('Shapiro-Wilk test statistics=%.3f, p=%.3f Sample looks Gaussian (fail to reject H0)' % (stat, p))
    else:
        normal_test_results[0] = 0
        print('Shapiro-Wilk test statistics=%.3f, p=%.3f Sample does not look Gaussian (reject H0)' % (stat, p))


    # D'Agostino's K**2 test 
    stat, p = normaltest(data)
    if p > alpha:
        print("D'Agostino's K**2 test statistics=%.3f, p=%.3f Sample looks Gaussian (fail to reject H0)" % (stat, p))
    else:
        normal_test_results[1] = 0
        print("D'Agostino's K**2 test statistics=%.3f, p=%.3f Sample does not look Gaussian (reject H0)" % (stat, p))

    # Anderson-Darling test
    result = anderson(data)
    if result[2][2] < result[1][2]:
        print('Anderson-Darling test statistic: %.3f, Sample looks Gaussian (fail to reject H0) sl: %.3f, cv: %.3f' % (result[0], result[2][2],result[1][2]))
    else:
        normal_test_results = 0
        print('Anderson-Darling test statistic: %.3f, Sample does not look Gaussian (reject H0) sl: %.3f, cv: %.3f' % (result[0], result[2][2],result[1][2]))
    # for i in range(len(result[1])):#critical_values)):
    #     sl, cv = result[2][i], result[1][i]
    #     if result[0] < result[1][i]:
    #         print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
    #     else:
    #         print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))

    is_normal = (np.mean(normal_test_results)>1)
    return is_normal

def run_stats(data1,data2,is_normal=False):
    data = data1-data2

    #tests of normality 
    #is_normal = test_normality(data)

    
    if is_normal: #parametric
        #print "\nParametric"
        cohens_d = (mean(data1) - mean(data2)) / (sqrt((stdev(data1) ** 2 + stdev(data2) ** 2) / 2)) #effect size power analysis
        t, p = ttest_1samp(data,0) #ttest
        print("Cohens d", cohens_d)
        print("ttest: t ", t, "p", p/2) #one-sided

    else: #nonparametric
        #print "\nNonparametric"
        #t,p = mannwhitneyu(data1,data2)
        #print "Mann-Whitney U: t ",t,"p", p
        t, p = wilcoxon(data1,data2)
        print("Wilcoxon: t ",t,"p", p)
        #print "Binomial sign test": np.sum((data1-data2)>0)/float(np.sum((data1-data2)!=0))
        p,data = resampling_statistics(data,0)
        print("p: ", p)
    return p
        
def bar_witherror_anddots(ax,data0,data1,x=[0,1],c0=[0/255.,0/255.,0/255.],c1=[0/255.,0/255.,0/255.]):
    
    n = np.size(data1)
    ax.scatter(np.zeros(n)+x[0],data0,s=100,facecolor='gray',edgecolor='None',alpha=.25,clip_on=False)
    ax.scatter(np.zeros(n)+x[1],data1,s=100,facecolor='gray',edgecolor='None',alpha=.25,clip_on=False)
    ax.plot(x,[data0,data1],color=[.5,.5,.5],alpha=.15)
    #ax.plot(x,[np.mean(data0),np.mean(data1)],linewidth=3,color='k')
    ax.bar(x[0],np.mean(data0),.5,linewidth=4,facecolor='None',edgecolor=c0)
    ax.bar(x[1],np.mean(data1),.5,linewidth=4,facecolor='None',edgecolor=c1)
    ax.errorbar(x[0],np.mean(data0),np.std(data0)/np.sqrt(n),color=c0,
                   lw=4,capsize=10,capthick=4,zorder=20)
    ax.errorbar(x[1],np.mean(data1),np.std(data1)/np.sqrt(n),color=c1,
                   lw=4,capsize=10,capthick=4,zorder=20)


def run_stats_onetail(data1,data2=None):
    '''
    run_stats_onetail runs parametric and nonparametric statistics for the paired comparison data1 > data2
    
    - data2 can either be a vector of the same size as data1, for a paired comparison
            or a single value, which means that it is chance
    '''

    if data2 is None:
        data = data1
        t, p = ttest_1samp(data1,0) #ttest
    elif np.size(data1)==np.size(data2):
        data = data1-data2
        t, p = ttest_rel(data1,data2) #ttest
    else:
        data = data1-data2
        t, p = ttest_1samp(data1,data2) #ttest

    #parametric - note: I didn't report parametric tests
    d = np.round(np.mean(data)/np.std(data),decimals=2)
    print("Parametric: ttest: t ", np.round(t,decimals=2), "p", '{:0.2e}'.format(p/2), "Cohens dz: ", d) #one-sided

    #nonparametric
    p,m,sd = resampling_statistics(data,0)
    if np.round(p,decimals=3)<0.001:
        print("Nonparametric p < 0.001")
    else:
        print("Nonparametric p: = {:.3f}".format(p))

    return p, d

def run_stats_twotail(data1,data2=None):
    '''
    run_stats_twotail runs parametric and nonparametric statistics for the paired comparison data1 != data2
    
    - data2 can either be a vector of the same size as data1, for a paired comparison
            or a single value, which means that it is chance
    '''
    
    if data2 is None:
        data = data1
    elif np.size(data1)==np.size(data2):
        data = data1-data2
    else:
        data = data1-data2
    
    #parametric - note: I didn't report parametric tests
    t, p = ttest_rel(data1,data2) #ttest
    d = np.round(np.mean(data)/np.std(data),decimals=2)
    print("Parametric: ttest: t ", np.round(t,decimals=2), "p", '{:0.2e}'.format(p), "Cohens dz: ", d) #one-sided

    #nonparametric
    p,m,sd = resampling_statistics(data,0)
    if p>.5:
        p = 1-p
    p = p*2
    if np.round(p,decimals=3)<0.001:
        print("Nonparametric p < 0.001")
    else:
        print("Nonparametric p: = {:.3f}".format(p))
    
    return p, d

def run_stats_btwnsubj(data1,data2):
    '''
    run_stats_onetail runs parametric and nonparametric statistics for the unpaired comparison data1 > data2
    
    - data1 and data2 do not need to be the same size
    '''
    
    #parametric
    t, p = ttest_ind(data1,data2) #ttest
    print("Parametric: ttest: t ", np.round(t,decimals=2), "p", '{:0.2e}'.format(p/nd))

    #nonparametric
    p,m1,m2 = resampling_statistics_btwnsubj(data1,data2)
    if np.round(p,decimals=3)<0.001:
        print("Nonparametric p < 0.001")
    else:
        print("Nonparametric p: = {:.3f}".format(p))

    return p

def bin_classify(X,y,clf,t=np.arange(-300,1501),t0 = np.arange(-300,1501,10),tw=100,nits=100,nbins=3,iclass=1,shuff=True):
    
    acc_bin = np.empty((nits,np.size(t0)))
    acc_bin_shuff = np.empty((nits,np.size(t0)))
    
#     if iclass==1:
#         clf = LogisticRegression(solver='liblinear',multi_class='ovr')
#     elif iclass==2:
#         clf = SVC(gamma=2, C=1)
#     elif iclass==3:
#         clf = SVC(kernel="linear", C=0.025)
#     elif iclass==4:
#         clf = KNeighborsClassifier(3)
#     elif iclass==5:
#         clf = GaussianProcessClassifier(1.0 * RBF(1.0))
#     elif iclass==6:
#         clf = DecisionTreeClassifier(max_depth=5)
#     elif iclass==7:
#         clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
#     elif iclass==8:
#         clf = MLPClassifier(alpha=1, max_iter=10000)
#     elif iclass==9:
#         clf = AdaBoostClassifier()
#     elif iclass==10:
#         clf = GaussianNB()
#     elif iclass==11:
#         clf = QuadraticDiscriminantAnalysis()

    cv = StratifiedShuffleSplit(n_splits=nits)
    
    nlabels=np.size(np.unique(y))
    nmin = np.min([np.sum(y==i) for i in np.arange(0,nlabels)])
    nperbin = int(np.floor(nmin/nbins))
    nmin = nperbin*nbins
    for it_counter in range(nits):

        y_bin = np.sort(np.tile((np.unique(y)),nbins))
        
        #balance sizes 
        X_bin = np.zeros((np.size(y_bin),np.shape(X)[1],np.shape(X)[2]))

        counter = 0
        for i,ilabel in enumerate(np.unique(y)):
            temp=np.where(y==ilabel)[0]
            np.random.shuffle(temp)            
            temp_binned_labels = np.reshape(temp[:nmin],(nbins,nperbin))
            for ibin in range(nbins):
                X_bin[counter] = np.mean(X[temp_binned_labels[ibin]],axis=0)
                counter=counter+1
        
        if nbins==2:
            itrain = np.sort(np.arange(0,nlabels*nbins,nbins))
            itest = np.arange(1,nlabels*nbins,nbins)            
        elif nbins==3:
            itrain = np.sort(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)))
            itest = np.arange(2,nlabels*nbins,nbins)
        elif nbins==4:
            itrain = np.sort(np.append(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)),np.arange(2,nlabels*nbins,nbins)))
            itest = np.arange(3,nlabels*nbins,nbins)
        elif nbins==6: 
            itrain = np.sort(np.append(np.append(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)),np.arange(2,nlabels*nbins,nbins)),np.arange(3,nlabels*nbins,nbins)))
            itest = np.sort(np.append(np.arange(4,nlabels*nbins,nbins),np.arange(5,nlabels*nbins,nbins)))
        else:
            print("only coded for nbins 2, 3, 4 or 6")
        
        for n,tstart in enumerate(t0):
            #find the indices of the appropriate time window
            idx = np.where(np.logical_and(t>=tstart,t<(tstart+tw)))[0]
            X_train = np.mean(X_bin[itrain][:,:,idx],axis=2)
            X_test = np.mean(X_bin[itest][:,:,idx],axis=2)

            y_train = y_bin[itrain]
            y_test = y_bin[itest]
            if shuff:
                y_train_shuff = y_bin[itrain][np.random.permutation(np.size(y_train))]
                y_test_shuff = y_bin[itest][np.random.permutation(np.size(y_test))]

            #train model
            clf.fit(X_train,y_train)
            acc_bin[it_counter,n] = clf.score(X_test,y_test)
            if shuff:
                clf.fit(X_train,y_train_shuff)
                acc_bin_shuff[it_counter,n] = clf.score(X_test,y_test_shuff)    

    return acc_bin, acc_bin_shuff

def singletrial_classify(X,y,clf,t=np.arange(-300,1501),t0 = np.arange(-300,1501,10),tw=100,nits=100,nbins=3,shuff=True):
    #NEED TO TEST THAT THIS FUNCTION WORKS! 
    acc = np.empty((nits,np.size(t0)))
    acc_shuff = np.empty((nits,np.size(t0)))
    
#     if iclass==1:
#         clf = LogisticRegression(solver='liblinear',multi_class='ovr')
#     elif iclass==2:
#         clf = SVC(gamma=2, C=1)
#     elif iclass==3:
#         clf = SVC(kernel="linear", C=0.025)
#     elif iclass==4:
#         clf = KNeighborsClassifier(3)
#     elif iclass==5:
#         clf = GaussianProcessClassifier(1.0 * RBF(1.0))
#     elif iclass==6:
#         clf = DecisionTreeClassifier(max_depth=5)
#     elif iclass==7:
#         clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
#     elif iclass==8:
#         clf = MLPClassifier(alpha=1, max_iter=10000)
#     elif iclass==9:
#         clf = AdaBoostClassifier()
#     elif iclass==10:
#         clf = GaussianNB()
#     elif iclass==11:
#         clf = QuadraticDiscriminantAnalysis()

    cv = StratifiedShuffleSplit(n_splits=nits)
    
    for itrain,itest in cv.split(X=X, y=y):
        
        for n,tstart in enumerate(t0):
            #find the indices of the appropriate time window
            idx = np.where(np.logical_and(t>=tstart,t<(tstart+tw)))[0]
            X_train = np.mean(X[itrain][:,:,idx],axis=2)
            X_test = np.mean(X[itest][:,:,idx],axis=2)

            y_train = y[itrain]
            y_test = y[itest]
            if shuff:
                y_train_shuff = y[itrain][np.random.permutation(np.size(y_train))]
                y_test_shuff = y[itest][np.random.permutation(np.size(y_test))]

            #train model
            clf.fit(X_train,y_train)
            acc[it_counter,n] = clf.score(X_test,y_test)
            if shuff:
                clf.fit(X_train,y_train_shuff)
                acc_shuff[it_counter,n] = clf.score(X_test,y_test_shuff)    
        
        it_counter = it_counter+1
        
    return acc, acc_shuff

def plot_classify(acc,acc_shuff,t,c1=[11/255.,0/255.,255/255.],c2='gray',nsubplots=2):

    nsubj = np.shape(acc)[0]
    
    if nsubplots==2:
        fig,ax = plt.subplots(1,2,figsize=(14,6))
        ax0 = ax[0]
        ax1 = ax[1]
    else:
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        ax0 = ax
    #plot classification accuracy
    y=np.mean(np.mean(acc,axis=1)*100,axis=0)
    yerr=np.std(np.mean(acc,axis=1)*100,axis=0)/np.sqrt(nsubj)
    ax0.fill_between(t,y-yerr,y+yerr,color=c1,edgecolor='None',alpha=.25,lw=0)
    ax0.plot(t,y,lw=4,color=c1)
    
    #plot shuffled null
    y=np.mean(np.mean(acc_shuff*100,axis=0),axis=0)
    yerr=np.std(np.mean(acc_shuff*100,axis=1),axis=0)/np.sqrt(nsubj)
    c=[100/255.,143/255.,255/255.]
    ax0.fill_between(t,y-yerr,y+yerr,color=c2,edgecolor='None',alpha=.25,lw=0)
    ax0.plot(t,y,color=c2,lw=4)
    
    #plot statistics
    pval = np.empty(np.size(t))
    for it in range(np.size(t)):
       _,pval[it] = (ttest_rel(np.mean(acc,axis=1)[:,it],np.mean(acc_shuff,axis=1)[:,it]))
    if np.mean(acc_shuff)<.20:
        ystat=29
    elif np.mean(acc_shuff)<.40:
        ystat=49
    else:
        ystat=69
    ax0.scatter(t[pval<.05],np.ones(np.sum(pval<.05))*ystat,color='k',clip_on=False,zorder=20)
    ax0.plot([0,0],[0,ystat+1],'--',color='k')
    
    #prettify plot
    if np.mean(acc_shuff)<.20:
        ax0.plot(t,np.mean(1/8.)*100*np.ones(np.size(t)),'--',color='k')
        prettify_plot(ax0,
                    xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                    ylim=[0,30],yt=[0,10,20,30],ytl=[0,10,20,30],yl='Classification accuracy (%)')
    elif np.mean(acc_shuff)<.40:
        ax0.plot(t,np.mean(1/4.)*100*np.ones(np.size(t)),'--',color='k')
        prettify_plot(ax0,
                    xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                    ylim=[20,50],yt=[20,30,40,50],ytl=[20,30,40,50],yl='Classification accuracy (%)')    
    else:
        ax0.plot(t,np.mean(1/2.)*100*np.ones(np.size(t)),'--',color='k')
        prettify_plot(ax0,
                xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                ylim=[40,70],yt=[40,50,60,70],ytl=[40,50,60,70],yl='Classification accuracy difference (%)')    
    
    #plot individual subjects
    if nsubplots==2:
        ax1.plot(t,np.zeros(np.size(t)),'--',color='k')
        ax1.plot(t,(np.mean(acc*100,axis=1)-np.mean(acc_shuff*100,axis=1)).T,color='k',alpha=.25)
        ax1.plot(t,np.mean((np.mean(acc*100,axis=1)-np.mean(acc_shuff*100,axis=1)),axis=0),color='k',lw=4)
        ax1.plot([0,0],[-20,40],'--',color='k')
        ax1.scatter(t[pval<.05],np.ones(np.sum(pval<.05))*(38+1/3.),color='k',clip_on=False)

        #prettify plot
        if np.mean(acc_shuff)<.20:
            prettify_plot(ax[1],
                    xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                    ylim=[-10,40],yt=[-40,-20,0,20,40],ytl=[-40,-20,0,20,40],yl='Classification accuracy difference (%)')
        elif np.mean(acc_shuff)<.40:
            prettify_plot(ax[1],
                    xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                    ylim=[-20,40],yt=[-40,-20,0,20,40],ytl=[-40,-20,0,20,40],yl='Classification accuracy difference (%)')    
        else:
            prettify_plot(ax[1],
                    xlim=[-300,1500],xt=[-300,0,300,600,900,1200,1500],xtl=[-300,0,300,600,900,1200,1500],
                    ylim=[-20,40],yt=[-40,-20,0,20,40],ytl=[-40,-20,0,20,40],yl='Classification accuracy difference (%)')    

    return fig


def run_nonparam_stats_onetail(data1,data2=None):
    '''
    run_stats_onetail runs parametric and nonparametric statistics for the paired comparison data1 > data2
    
    - data2 can either be a vector of the same size as data1, for a paired comparison
            or a single value, which means that it is chance
    '''

    if data2 is None:
        data = data1
        #t, p = ttest_1samp(data1,0) #ttest
    elif np.size(data1)==np.size(data2):
        data = data1-data2
        #t, p = ttest_rel(data1,data2) #ttest
    else:
        data = data1-data2
        #t, p = ttest_1samp(data1,data2) #ttest

    #parametric - note: I didn't report parametric tests
#    d = np.round(np.mean(data)/np.std(data),decimals=2)
#    print("Parametric: ttest: t ", np.round(t,decimals=2), "p", '{:0.2e}'.format(p/2), "Cohens dz: ", d) #one-sided

    #nonparametric
    p,m,sd = resampling_statistics(data,0)
#     if np.round(p,decimals=3)<0.001:
#         print("Nonparametric p < 0.001")
#     else:
#         print("Nonparametric p: = {:.3f}".format(p))

    return p


def bin_classify_probcoef(X,y,clf,t=np.arange(-300,1501),t0 = np.arange(-300,1501,10),tw=100,nits=100,nbins=3,iclass=1,shuff=True,scaler=True,trainprop=.5):
    
    acc_bin = np.empty((nits,np.size(t0)))
    coef_bin = np.empty((nits,np.size(t0),np.size(np.unique(y)),np.shape(X)[1]))
    proba_bin = np.empty((nits,np.size(t0),np.size(np.unique(y)),np.size(np.unique(y))))
    acc_bin_shuff = np.empty((nits,np.size(t0)))
    
    cv = StratifiedShuffleSplit(n_splits=nits)
    scaler = StandardScaler()
    
    nlabels=np.size(np.unique(y))
    nmin = np.min([np.sum(y==i) for i in np.arange(0,nlabels)])
    nperbin = int(np.floor(nmin/nbins))
    nmin = nperbin*nbins
    for it_counter in range(nits):

        y_bin = np.sort(np.tile((np.unique(y)),nbins))
        
        #balance sizes 
        X_bin = np.zeros((np.size(y_bin),np.shape(X)[1],np.shape(X)[2]))

        counter = 0
        for i,ilabel in enumerate(np.unique(y)):
            temp=np.where(y==ilabel)[0]
            np.random.shuffle(temp)
            temp = temp[:nmin] #downsample
            if nbins != 2:
                temp_binned_labels = np.reshape(temp,(nbins,nperbin))
            else:
                if trainprop==.5:
                    temp_binned_labels = np.reshape(temp,(nbins,nperbin))
                else:
                    temp_binned_labels = {}
                    i = int(np.floor(trainprop*nmin))
                    temp_binned_labels[0] = temp[:i]
                    temp_binned_labels[1] = temp[i:]
            for ibin in range(nbins):
                X_bin[counter] = np.mean(X[temp_binned_labels[ibin]],axis=0)
                counter=counter+1
        
        if nbins==2:
            itrain = np.sort(np.arange(0,nlabels*nbins,nbins))
            itest = np.arange(1,nlabels*nbins,nbins)            
        elif nbins==3:
            itrain = np.sort(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)))
            itest = np.arange(2,nlabels*nbins,nbins)
        elif nbins==4:
            itrain = np.sort(np.append(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)),np.arange(2,nlabels*nbins,nbins)))
            itest = np.arange(3,nlabels*nbins,nbins)
        elif nbins==6: 
            itrain = np.sort(np.append(np.append(np.append(np.arange(0,nlabels*nbins,nbins),np.arange(1,nlabels*nbins,nbins)),np.arange(2,nlabels*nbins,nbins)),np.arange(3,nlabels*nbins,nbins)))
            itest = np.sort(np.append(np.arange(4,nlabels*nbins,nbins),np.arange(5,nlabels*nbins,nbins)))
        else:
            print("only coded for nbins 2, 3, 4 or 6")
        
        for n,tstart in enumerate(t0):
            #find the indices of the appropriate time window
            idx = np.where(np.logical_and(t>=tstart,t<(tstart+tw)))[0]
            X_train = np.mean(X_bin[itrain][:,:,idx],axis=2)
            X_test = np.mean(X_bin[itest][:,:,idx],axis=2)
            
            #standard scaler
            scaler.fit(X_train)
            scaler.transform(X_train)
            scaler.transform(X_test)
            
            y_train = y_bin[itrain]
            y_test = y_bin[itest]
            if shuff:
                y_train_shuff = y_bin[itrain][np.random.permutation(np.size(y_train))]
                y_test_shuff = y_bin[itest][np.random.permutation(np.size(y_test))]

            #train model
            clf.fit(X_train,y_train)
            
            #test model
            acc_bin[it_counter,n] = clf.score(X_test,y_test)
            
            #coef
            #print(clf.coef_)
            #print(np.shape(clf.coef_))
            coef_bin[it_counter,n] = clf.coef_
            
            #proba
            proba_bin[it_counter,n] = clf.predict_proba(X_test)
            if shuff:
                clf.fit(X_train,y_train_shuff)
                acc_bin_shuff[it_counter,n] = clf.score(X_test,y_test_shuff)    

    return acc_bin, coef_bin, proba_bin, acc_bin_shuff
