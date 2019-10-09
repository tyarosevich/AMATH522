import numpy as np
from datetime import datetime
import random
import os.path
import time
import json
import gzip

version = 'final_forRepository'
subDirString = 'K8init500Tratio1'
saveToFile = True
saveGridsInEachGen = False # normally set to False. If True will save the Bact lattice in each gen

stopFrac = np.nan # set to np.nan for fig3; and to 0.05 for fig4

N = 100       # lattice size: N X N
gen = 10000   # number of generations in simulation
nei = 1       # neighbourhood radius. 1 -> 8 neighbors. paper is entirely with nei=1

# params for initializations of latices
baseLineC = 0
initialNumCoopBact = int(0.05*N*N) # set int(0.05*N*N) for fig3; set to int(4) for fig4
coopBactStartInPatch = False # if False distribute randomally; False for fig3, True for fig4 

# initialNumA refers to number of individuals who choose neighbors to interact with in a non-random fashion
# initialNumA must be 0 since currently code only implements random choice of neighbor for interaction
initialNumA = 0 # implemented only for 0; don't set to a different value

# list of parameter values to choose from
KValues = [8]
bOverCValues = np.concatenate((np.arange(1, 5.5, 0.5), [10, 100]))
cIncValues = [0.05]
infectionValues = np.arange(0, 0.525, 0.05) 
infectionRatioValues = [1]
verticalNoise = 0 # VT = 1-verticalNoise; when 0 there is no noise

if saveGridsInEachGen: # usually used for a single paramset. Below steps over usual parameters
    KValues = [1]
    bOverCValues = [3]
    cIncValues = [0.05]
    infectionValues = [0.2]
    infectionRatioValues = [1]

windowForSmoothing = 200

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyEncoder, self).default(obj)

def neighbours(mat, x, y, nx=1, ny=1):
    ''' returns a slice from a matrix (mat) with radius (nx,ny) around coordinate x,y
    nx,ny get difult value 1''' 
    return mat[max(x-nx, 0) : min(x+nx+1, mat.shape[0]), \
                      max(y-ny, 0) : min(y+ny+1, mat.shape[1])]

def find_ind_rand(m,mat):
    '''returns at random one of the indices in matrix (mat) of a specific value (m)'''
    found = np.argwhere(np.isclose(mat,m))
    r = np.random.randint(0,found.shape[0]) 
    return found[r,0],found[r,1]
    
def ind2sub(numCols,ind):
    row = (ind // numCols)
    col = (ind % numCols) 
    return (row, col)

def interactions(laticeC,laticeAa,laticeBact,params):
    '''calculates the fitness of each individual according to the chosen interactors out of its neighbours
    the function doesn't return a value, but rather changes the fitness matrix (laticefi) instead.
    The function also changes the bacteria matrix due to interactions'''
    N = params['N']
    nei = params['nei']
    laticeIds = np.arange(N*N).reshape(N,N)
    laticefi = np.zeros((N,N), dtype=float)
    laticeNumInteractions = np.zeros((N,N))   
    for it in range(params['K']):
        allSubs = [ (i,j) for i in range(N) for j in range(N) ] 
        random.shuffle(allSubs)
        for mySub in allSubs:
            i = mySub[0]
            j = mySub[1]
            myId = i*N+j 
            neibIds = neighbours(laticeIds,i,j,nei,nei).flatten().tolist()
            neibIds.remove(myId) # remove self from list of neighbor ids
            neibId = random.choice(neibIds)
            neibSub = ind2sub(N,neibId)
            myfi, neibfi = individualInteraction(mySub,neibSub,laticeC,laticeBact, params)            
            laticeNumInteractions[i,j] += 1
            laticeNumInteractions[neibSub] += 1              
            laticefi[i,j] += myfi
            laticefi[neibSub] += neibfi
    laticefi /= laticeNumInteractions
    return laticefi

def individualInteraction(mySub,neibSub,laticeC,laticeBact,params):
    # returns fitness delta due to interaction and changes bacteria phentype for interactors
    cInc = params['cInc']   
    infection = params['infection']
    infectionRatio = params['infectionRatio']
    myBact = laticeBact[mySub]
    neibBact = laticeBact[neibSub]
    myC = laticeC[mySub] + myBact*cInc # if Bact=0 no addition to coop, if bact=1 adds cInc
    neibC = laticeC[neibSub] + neibBact*cInc
    myfi = neibC*bOverC - myC
    neibfi = myC*bOverC - neibC
    if neibBact != myBact:
        if neibBact == 1: # myBact is 0
            if np.random.rand()<infection: # cooperative bacteria infects by infection
                laticeBact[mySub] = 1
            if np.random.rand()<infection*infectionRatio: # non-cooperative bacteria infects by infection*infectionRatio
                laticeBact[neibSub] = 0
        else: # myBact is 1; neiBact is 0
            if np.random.rand()<infection*infectionRatio:
                laticeBact[mySub] = 0
            if np.random.rand()<infection:
                laticeBact[neibSub] = 1
    return myfi, neibfi

def calcNeibStats(laticeBact):
    # returns dict of number of 0-0/0-1/1-0/1-1 pairs of neighbors
    numBact0Bact0Neibs = 0   
    numBact0Bact1Neibs = 0 # symmetric
    numBact1Bact1Neibs = 0 
    
    rolledDown = np.roll(laticeBact,1,axis=0)
    b = laticeBact[1:,:] + rolledDown[1:,:]
    numBact0Bact0Neibs += np.sum(b==0)
    numBact0Bact1Neibs += np.sum(b==1)
    numBact1Bact1Neibs += np.sum(b==2)

    rolledRight = np.roll(laticeBact,1,axis=1)
    b = laticeBact[:,1:] + rolledRight[:,1:]
    numBact0Bact0Neibs += np.sum(b==0)
    numBact0Bact1Neibs += np.sum(b==1)
    numBact1Bact1Neibs += np.sum(b==2)
    
    rolledDownAndRight = np.roll(np.roll(laticeBact,1,axis=1),1,axis=0)
    b = laticeBact[1:,1:] + rolledDownAndRight[1:,1:]
    numBact0Bact0Neibs += np.sum(b==0)
    numBact0Bact1Neibs += np.sum(b==1)
    numBact1Bact1Neibs += np.sum(b==2)
    
    rolledDownAndLeft = np.roll(np.roll(laticeBact,-1,axis=1),1,axis=0)
    b = laticeBact[1:,0:-1] + rolledDownAndLeft[1:,0:-1]
    numBact0Bact0Neibs += np.sum(b==0)
    numBact0Bact1Neibs += np.sum(b==1)
    numBact1Bact1Neibs += np.sum(b==2)
    
    neibStats = {
            'numBact0Bact0Neibs': numBact0Bact0Neibs,   
            'numBact0Bact1Neibs': numBact0Bact1Neibs,   
            'numBact1Bact1Neibs': numBact1Bact1Neibs,   
            }
    return neibStats
                    
def simulation(laticeC,laticeAa,laticeBact,params):
    N = params['N']
    nei = params['nei']
    verticalNoise = params['verticalNoise']
    saveGridsInEachGen = params['saveGridsInEachGen']
    startTime = time.clock()
    fracBact = [] 
    fracBactSmoothed = []
    fracsPreInteractionPerGeneration = []
    fracsPostInteractionPerGeneration = []
    meanFitnessPreInteractionPerGeneration = []
    meanFitnessPostInteractionPerGeneration = []
    neibStatsPerGeneration = []
    
    if saveGridsInEachGen:
        writeEveryGen = 200 # collect data, and export to file every writeEveryGen generations
        laticeBactPreInteractionPerGen = np.zeros((writeEveryGen*N,N))-1 # keeping it a 2d-array for savetxt. N*N matrices one under the other
        laticeBactPostInteractionPerGen = np.zeros((writeEveryGen*N,N))-1 
        laticeFitnessPerGen = np.zeros((writeEveryGen*N,N))-1 
        fnPreInteractionBactBase = params['namePreInteractionBactFileFullPath']
        fnPostInteractionBactBase = params['namePostInteractionBactFileFullPath']
        fnFitnessBase = params['nameFitnessFileFullPath']
        part = 1
        
    for k in range(gen):   # loop over generations
        # check stop criteria
        if k>0:
            if np.isclose(fracBact[-1],0) or np.isclose(fracBact[-1],1):
                endTime = time.clock()
                print('Finished due to 0/1 criterion')
                break
            if ~np.isnan(params['stopFrac']) and fracBact[-1]>params['stopFrac']:
                endTime = time.clock()
                print('Finished due to stopFrac criterion')
                break
            if k%500 == 0:
                if abs(np.max(fracBactSmoothed[-200:])-np.min(fracBactSmoothed[-200:]))<0.01:
                    endTime = time.clock()
                    print('Finished due to constant smoothed fracBact criterion')     
                    break
            
        newC = np.zeros(laticeC.shape)
        newAa = np.zeros(laticeAa.shape)
        newBact = np.zeros(laticeBact.shape)
        
        # at this stage we have the bact status pre interaction ("at birth")
        fracBact.append(float(np.sum(laticeBact))/(N*N))
        fracBactSmoothed.append(np.mean(fracBact[max(0,k-params['windowForSmoothing']):k+1]))
        fracsPreInteractionPerGeneration.append({
            'bact0A0': float(np.sum((laticeBact==0)*(laticeAa==0))) / (N*N),
            'bact0A1': float(np.sum((laticeBact==0)*laticeAa)) / (N*N),
            'bact1A0': float(np.sum(laticeBact*(laticeAa==0))) / (N*N),
            'bact1A1': float(np.sum(laticeBact*laticeAa)) / (N*N),
        })
        neibStatsPerGeneration.append(calcNeibStats(laticeBact))

        # perform the interactions
        laticeBactPreInteractions = laticeBact.copy() # save a copy so that can calc mean fitness based on pre-interactions bact
        laticefi = interactions(laticeC,laticeAa,laticeBact,params)
        
        # calc the mean fitness using the pre interactions bacts
        meanFitnessPreInteractionPerGeneration.append({
            'bact0A0': float(np.sum((laticeBactPreInteractions==0)*(laticeAa==0)*laticefi)) / (N*N),
            'bact0A1': float(np.sum((laticeBactPreInteractions==0)*laticeAa*laticefi)) / (N*N),
            'bact1A0': float(np.sum(laticeBactPreInteractions*(laticeAa==0)*laticefi)) / (N*N),
            'bact1A1': float(np.sum(laticeBactPreInteractions*laticeAa*laticefi)) / (N*N),
        })
        
        # at this stage we have the bact status post interaction and their fitnesses after interactions
        fracsPostInteractionPerGeneration.append({
            'bact0A0': float(np.sum((laticeBact==0)*(laticeAa==0))) / (N*N),
            'bact0A1': float(np.sum((laticeBact==0)*laticeAa)) / (N*N),
            'bact1A0': float(np.sum(laticeBact*(laticeAa==0))) / (N*N),
            'bact1A1': float(np.sum(laticeBact*laticeAa)) / (N*N),
        })
        meanFitnessPostInteractionPerGeneration.append({
            'bact0A0': float(np.sum((laticeBact==0)*(laticeAa==0)*laticefi)) / (N*N),
            'bact0A1': float(np.sum((laticeBact==0)*laticeAa*laticefi)) / (N*N),
            'bact1A0': float(np.sum(laticeBact*(laticeAa==0)*laticefi)) / (N*N),
            'bact1A1': float(np.sum(laticeBact*laticeAa*laticefi)) / (N*N),
        })

        if saveGridsInEachGen:
            if k>0 and k%writeEveryGen==0:
                with open(fnPostInteractionBactBase+'part'+str(part),'wb') as f:
                    np.savetxt(f, laticeBactPostInteractionPerGen)
                laticeBactPostInteractionPerGen = np.zeros((writeEveryGen*N,N))-1
                with open(fnPreInteractionBactBase+'part'+str(part),'wb') as f:
                    np.savetxt(f, laticeBactPreInteractionPerGen)
                laticeBactPreInteractionPerGen = np.zeros((writeEveryGen*N,N))-1
                with open(fnFitnessBase+'part'+str(part),'wb') as f:
                    np.savetxt(f, laticeFitnessPerGen)
                laticeFitnessPerGen = np.zeros((writeEveryGen*N,N))-1
                part += 1
            laticeBactPreInteractionPerGen[(k%writeEveryGen)*N:((k%writeEveryGen)+1)*N,:] = laticeBactPreInteractions
            laticeBactPostInteractionPerGen[(k%writeEveryGen)*N:((k%writeEveryGen)+1)*N,:] = laticeBact
            laticeFitnessPerGen[(k%writeEveryGen)*N:((k%writeEveryGen)+1)*N,:] = laticefi
            
        for i in range(N):
            for j in range(N):
                neibC = neighbours(laticeC,i,j,nei,nei)
                neibfi = neighbours(laticefi,i,j,nei,nei)
                neibAa = neighbours(laticeAa,i,j,nei,nei)
                neibBact = neighbours(laticeBact,i,j,nei,nei)
                m = np.max(neibfi)  # find highest fitness value in neighbourhood
                sub = find_ind_rand(m,neibfi)
                newC[i,j] = neibC[sub]
                newAa[i,j] = neibAa[sub]
                if verticalNoise > 0:
                    if np.random.random() < verticalNoise:
                        if np.random.random() < neibBact.mean():
                            newBact[i,j] = 1
                        else:
                            newBact[i,j] = 0
                    else:
                        newBact[i,j] = neibBact[sub]
                else:
                    newBact[i,j] = neibBact[sub]
        laticeC = newC
        laticeAa = newAa
        laticeBact = newBact 
        
    endTime = time.clock()
    print('Successfully finished after %d generations. Running time in seconds: %.2f' %(k,endTime-startTime) )
    if saveGridsInEachGen and k%writeEveryGen!=0:
        with open(fnPreInteractionBactBase+'part'+str(part),'wb') as f:
            np.savetxt(f, laticeBactPreInteractionPerGen)
        with open(fnPostInteractionBactBase+'part'+str(part),'wb') as f:
            np.savetxt(f, laticeBactPostInteractionPerGen)
        with open(fnFitnessBase+'part'+str(part),'wb') as f:
            np.savetxt(f, laticeFitnessPerGen)

    return laticeC, laticefi, laticeAa, laticeBact, \
        fracsPreInteractionPerGeneration, fracsPostInteractionPerGeneration, \
        meanFitnessPreInteractionPerGeneration, meanFitnessPostInteractionPerGeneration, \
        neibStatsPerGeneration
                
# choose one parameter set randomly 
K = np.random.choice(KValues)
bOverC = np.random.choice(bOverCValues)
cInc = np.random.choice(cIncValues)
infection = np.random.choice(infectionValues)
infectionRatio = np.random.choice(infectionRatioValues)

while (infection==0 and infectionRatio!=infectionRatioValues[0]):
    infection = np.random.choice(infectionValues)
    infectionRatio = np.random.choice(infectionRatioValues)
    
''' initialize lattice of C''' # C is a continous number
laticeC = baseLineC*np.ones((N,N), dtype=float)

''' initialize lattice of A/a''' # Implemented only for random choice of neighbor (LaticeAa is all 0)
assert initialNumA==0
arr = np.array([1] * initialNumA + [0] * (N*N-initialNumA), dtype = int)
np.random.shuffle(arr)
laticeAa = arr.reshape(N,N)

''' initialize lattice of bacteria''' # 0/1, 1 is the cooperative bacteria
if coopBactStartInPatch:
    assert np.isclose(np.sqrt(initialNumCoopBact)%1,0) , 'initialNumCoopBact must be a square of an int'
    d = int(np.sqrt(initialNumCoopBact))
    t = int(np.floor((N-d)/2))
    laticeBact = np.zeros((N,N), dtype=int)
    laticeBact[t:t+d,t:t+d] = 1
else:
    arr = np.array([1] * initialNumCoopBact + [0] * (N*N-initialNumCoopBact), dtype = int)
    np.random.shuffle(arr)
    laticeBact = arr.reshape(N,N)   

### file names ###
now = datetime.now()
date = now.strftime('%Y-%b-%d_%H-%M-%S-%f')
nameOfile = 'nowak_spatial_modified_v' + str(version)+ '_out_{0}_.json.gz'.format(date)
if saveGridsInEachGen:
    namePreInteractionBactFile = 'nowak_spatial_modified_v' + str(version)+ '_preInteractionBact_{0}_'.format(date)
    namePostInteractionBactFile = 'nowak_spatial_modified_v' + str(version)+ '_postInteractionBact_{0}_'.format(date)
    nameFitnessFile = 'nowak_spatial_modified_v' + str(version)+ '_fitness_{0}_'.format(date)

### set params dict ###
params = {
            'N': N,
            'nei': nei,
            'gen': gen,
            'verticalNoise': verticalNoise,
            'initialNumA': initialNumA,
            'initialNumCoopBact': initialNumCoopBact,
            'coopBactStartInPatch': coopBactStartInPatch,
            'baseLineC': baseLineC,
            'K': K,
            'bOverC': bOverC,
            'cInc': cInc,
            'infection': infection,
            'infectionRatio': infectionRatio,
            'windowForSmoothing': windowForSmoothing,
            'stopFrac': stopFrac, 
            'saveGridsInEachGen': saveGridsInEachGen,
}
if saveToFile:   
    savePath = 'runs/v'+version+'/'+subDirString
    if not os.path.exists(savePath):
        os.makedirs(savePath)
    fn = os.path.join(savePath, nameOfile)

if saveGridsInEachGen:
    fn_PreInteractionbact = os.path.join(savePath, namePreInteractionBactFile)
    fn_PostInteractionbact = os.path.join(savePath, namePostInteractionBactFile)
    fn_fitness = os.path.join(savePath, nameFitnessFile)
    params['namePreInteractionBactFileFullPath'] = fn_PreInteractionbact
    params['namePostInteractionBactFileFullPath'] = fn_PostInteractionbact
    params['nameFitnessFileFullPath'] = fn_fitness
    
print(params)
##############################

finalLaticeC, finalLaticefi, finalLaticeAa, finalLaticeBact, \
        fracsPreInteractionPerGeneration, fracsPostInteractionPerGeneration, \
        meanFitnessPreInteractionPerGeneration, meanFitnessPostInteractionPerGeneration, \
        neibStatsPerGeneration = simulation(laticeC,laticeAa,laticeBact,params)

if saveToFile:   
    jsonResults = params.copy()
    jsonResults['fracsPreInteractionPerGeneration'] = fracsPreInteractionPerGeneration
    jsonResults['fracsPostInteractionPerGeneration'] = fracsPostInteractionPerGeneration
    jsonResults['meanFitnessPreInteractionPerGeneration'] = meanFitnessPreInteractionPerGeneration
    jsonResults['meanFitnessPostInteractionPerGeneration'] = meanFitnessPostInteractionPerGeneration
    jsonResults['neibStatsPerGeneration'] = neibStatsPerGeneration
    openfile = gzip.open
    with openfile(fn, 'wt') as f:
        json.dump(jsonResults, f, sort_keys = True, indent = 4, separators = (',', ':'), cls=NumpyEncoder)
        