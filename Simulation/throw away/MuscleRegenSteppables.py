
from cc3d.core.PySteppables import *
#import pandas as pd
import numpy as np
import random
import math

#VARIABLES#

# Simulation Variables
# Dimensions need to be updated for every new histology image
xdim = 321
ydim = 417
timeconverter = 60/41 # 1 mcs = 41 minutes, conversion factor from hours to mcs
recruitMultiplier = 2 # Multiplier of number recruited cells to check locations
InitialNecroticPercentage = 0.3 # Percent of fiber that are initialized as necrotic

# Neutrophil Variables 
NeutrophilPhagocytoseLagTime = 10  # Lag time after arrival before neutrophils phagocytose
neutroRecruitmentFreq = 1          # Frequency to recruit neutrophils in steps
neutroRecruitmentProportion = 0.15 # Number of neutrophils generated at every recruitment as proportion of necrotic cells
MCPsecretionVal = 200              # MCP secreted by neutrophil per step during and after "eating necrosis"
MMPsecretionVal = 200              # MMP secreted by neutrophil per step during and after "eating necrosis"
#neuToInitialFiber = 1/450         # Used if including resident neurophils
TimeNeuStart = math.ceil(2 * timeconverter)       # Start time of neutrophil recruitment in steps
neutroStopRecruit = math.ceil(24 * timeconverter) # Stop time of neutrophil recruitment in steps
neutroDeathTime = math.ceil(12.5 * timeconverter) # Time after birth for neutrophil to die
neutrophilHGFUptakeMax = 20        # Max uptake of HGF by neutrophils
neutrophilHGFUptakeRel = 0.3       # Relative uptake of HGF by neutrophils

# Macrophage Variables 
macToInitialFiber = 1/3.7             # Used if including resident macrophages, num of macrophages is proportion of original fibers
macroRecruitmentFreq = 1              # Start time of macrophage recruitment in steps
origMacroRecruitmentProportion = 0.8  # Recruitment of macrophages is proportional to mean MCP
MacrophageDeathTime = 50              # Added to randomly generate number to get death time
MacrophageDeathPoisson = 70           # Used in poisson distribution to get randomness in death time
MacrophageRandomDeath = 75            # initially, 1 in 75 chance of dying, decrease by 1 every timestep

TimeMacroStart = math.ceil(24 * timeconverter)  # Start time for macrophage recruitment in steps
MacrophagePhagocytoseLagTime = math.ceil((20/60) * timeconverter) # Lag after recruitment to start phagocytosis
macroStopRecruit = math.ceil(4*24 * timeconverter) # Stop time for macrophage recruitment in steps
macrophageHGFUptakeMax = 20           # Max uptake of HGF by macrophages
macrophageHGFUptakeRel = 0.3          # Relative uptake of HGF by macrophages
macrophageMCPUptakeMax = 20           # Max uptake of MCP by macrophages
macrophageMCPUptakeMRel = 0.3         # Relative uptake of MCP by macrophages
macrophageE2UptakeMax = 20/100.0      # Max uptake of E2 by macrophages
macrophageE2UptakeRel = 0.3/100.0     # Relative uptake of E2 by macrophages

# SSC Variables 
sscActivationThreshold = 1.5    # Level of HGF to activate ssc
sscToInitialFiber = 1/256       # Starting ssc are propotional to number of fibers
recruitmentProportionSSC = 0.15 # Number of ssc recruited every time step is proportional to sum of hgf and mmp means
moveCountThreshold = 20         # Number of moves without seeing low collagen ecm required to deactivate
sscApoptosisCutoff = 30         # number of moves without low collagen ecm before triggering apoptosis
sscActivationTime = math.ceil(1.5 * timeconverter)     # Lag time for activation
sscDivisionTime = math.ceil(10*timeconverter)          # Lag time for division
sscDifferentiationTime = 26     # Lag time for differentiation to fiber


# Necrosis Variables 
necBreakdownLag = 50            # length of time after healing that neutrophils and macrophages secrete cytokines
lowCollagenCutoff = 0.9         # cutoff value between normal and low collagen ECM
HGFsecretionVal = 50            # amount of hgf secreted by necrotic cells per step # ADD TO OPTIMIZATION CODE

# Estrogen Variables 
estrogenOn = False # True to use estrogen function, False to exclude 
systemicE2 = 0.0340 # based on starting data in proestrus
macAlpha = 1
MCPAlpha = 1
NeutroRecruitementProAlpha = 0.05
MacroRecruitementProAlpha = 17
NecBreakdownLagAlpha = 1
sscApoptosisAlpha = 1
E2CycleIndex = 0
# EstrusCycle2 data from McLean et al., 2012 (time step for every 6hrs = 8.78 mcs)
cycleData = [0.0334, 0.0476, 0.0647, 0.0612, 0.0560, 0.0433, 0.0273, 0.0185, 0.0153, 0.0193, 0.0281, 0.0370, 0.0300, 0.0225, 0.0170, 0.0226, 0.0332]
sscE2UptakeMax = 20/100.0
sscE2UptakeRel = 0.3/100.0

# Lymphatic Variables 
lympForce = 100 #  Force pulling ssc, neutrophils, and macrophages into lymph

class ImportInitializerSteppable(SteppableBasePy): # Only run to segment histology image 
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        fileName = 'C:\\Users\\megan\\OneDrive - University of Virginia\\M3\\EEE\\CC3D\\Hackathon\\TestImport\\healthy_histo2.csv'
        histoData = pd.read_csv(fileName, header=None)
        upscale = True
        if not upscale:
            for i in range(len(histoData)):
                x = histoData.iloc[i,0]
                y = histoData.iloc[i,1]
                tissueType = histoData.iloc[i,2]
                
                # if tissueType == 1: # ?
                    # cell_ecm = self.new_cell(self.WALL)
                    # self.cellField[x.tolist(),y.tolist(),0] = cell_ecm
                    
                # if tissueType == 3: # Fiber
                    # cell_fib = self.new_cell(self.FIBER)
                    # self.cellField[x.tolist(),y.tolist(),0] = cell_fib

                if tissueType == 2: # ECM 
                    cell_ecm = self.new_cell(self.WALL)
                    self.cellField[x.tolist(),y.tolist(),0] = cell_ecm
                
        if upscale:
            for i in range(len(histoData)):
                x = histoData.iloc[i,0].tolist()
                y = histoData.iloc[i,1].tolist()
                tissueType = histoData.iloc[i,2]
                
                xAdjusted = (x-1)*3 +1
                yAdjusted = (y-1)*3 +1
                
                if tissueType == 1: # ?
                    cell_ecm = self.new_cell(self.WALL)
                    for i in range(xAdjusted, xAdjusted+3):
                        for j in range(yAdjusted, yAdjusted+3):
                            self.cellField[i,j,1] = cell_ecm
                    
                # if tissueType == 3: # Fiber
                    # cell_fib = self.new_cell(self.FIBER)
                    # for i in range(xAdjusted, xAdjusted+2):
                        # for j in range(yAdjusted, yAdjusted+2):
                            # self.cellField[i,j,0] = cell_fib

                # if tissueType == 2: # ECM 
                    # cell_fib = self.new_cell(self.ECM)
                    # for i in range(xAdjusted, xAdjusted+3):
                        # for j in range(yAdjusted, yAdjusted+3):
                            # self.cellField[i,j,0] = cell_fib
        ## each square to become 3 by 3
    
    def step(self,mcs):
        return
        # if mcs == 1:
            # foundSpace = False
            
            # while not foundSpace:
                # xrand = np.random.uniform(0, self.dim.x) 
                # yrand = np.random.uniform(0, self.dim.y)
               
                # cell = self.cellField[xrand,yrand,0]
                # if not cell:
                    # cell_fib = self.new_cell(self.FIBER)
                    # self.cellField[xrand, yrand,0] = cell_fib
                    # cell_fib.targetVolume = 9500
                    # cell_fib.lambdaVolume =1000
                    # foundSpace = True
                    
        
        # allFilled = True
        # for cell in self.cell_list_by_type(self.FIBER):
            # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                # if not neighbor:
                    # allFilled = False
                    
            
        # if allFilled:
            # foundSpace = False
            
            # while not foundSpace:
                # xrand = np.random.uniform(0, self.dim.x) 
                # yrand = np.random.uniform(0, self.dim.y)
               
                # cell = self.cellField[xrand,yrand,0]
                # if not cell:
                    # cell_fib = self.new_cell(self.FIBER)
                    # self.cellField[xrand, yrand,0] = cell_fib
                    # cell_fib.targetVolume = 9500
                    # cell_fib.lambdaVolume =1000
                    # foundSpace = True  

class MitosisSteppable(MitosisSteppableBase): # Only run to divide up fibers from orginial fiber segmentation image
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
    
    def step(self, mcs):
        if mcs == 1:
            i = 0
            while i < 6: #randomly divide fiber fasicles into smaller regions for necrosis function
                cells_to_divide=[]
                for cell in self.cell_list_by_type(self.FIBER):
                    if cell.volume > 30:
                        cells_to_divide.append(cell)

                for cell in cells_to_divide:
                    self.divide_cell_random_orientation(cell)
                    # Other valid options
                    # self.divide_cell_orientation_vector_based(cell,1,1,0)
                    # self.divide_cell_along_major_axis(cell)
                    # self.divide_cell_along_minor_axis(cell)
                i = i+1
                # cells_to_divide=[]
                # for cell in self.cell_list:
                    # if cell.volume>50000:
                        # cells_to_divide.append(cell)

                # for cell in cells_to_divide:

                    # self.divide_cell_random_orientation(cell)
                    # # Other valid options
                    # # self.divide_cell_orientation_vector_based(cell,1,1,0)
                    # # self.divide_cell_along_major_axis(cell)
                    # # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0 
        self.parent_cell.lambdaVolume/= 2.0
        self.parent_cell.fluctAmpl = 50
        self.clone_parent_2_child()        

        # for more control of what gets copied from parent to child use cloneAttributes function
        self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[])
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2])
        
        # if self.parent_cell.type==1:
            # self.child_cell.type=2
        # else:
            # self.child_cell.type=1


class ConstraintInitializerSteppable(SteppableBasePy): # Initialize properties of ecm and fiber imported with piff
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        # Initialize volume target and lambda for ecm and fiber imported with piff
        for cell in self.cell_list_by_type(self.FIBER):
            cell.targetVolume = cell.volume     # Set target volume equal to volume defined by piff
            cell.lambdaVolume = 10000           # Has to be large to have stiff fibers
        
        for cell in self.cell_list_by_type(self.ECM):
            cell.targetVolume = cell.volume + 1 # Set target volume to slighly above that defined by piff
            cell.lambdaVolume = 10000           # Has to be large to have stiff ecm
            cell.dict['collagen'] = 1.5         # Define collagen level above low collagen threshold
                                  
class NecrosisSteppable(SteppableBasePy): # Initialize muscle injury
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # calculate total volume in pixels of fibers before necrosis happens. loop through all fiber cells, add their volume in pixels to FiberVolume
        # size of cell will be 3x3x1
        return

    def step(self, mcs):
        # Gets random points and gets area around it to place and spread necrosis
        NumNecrosis = 100
        NecrosisRadius = 0.03 * self.dim.x
        
        # Generate necrosis at step 0
        if mcs == 0: 
            NecroticVolume = 0
            FiberVolume = 0
            for cell in self.cell_list_by_type(self.FIBER):
                FiberVolume+= cell.volume
            NecroticPercentage = NecroticVolume / FiberVolume
            
            # Generate array of random coordinates
            xarray = np.random.uniform(0,self.dim.x,NumNecrosis)
            yarray = np.random.uniform(0,self.dim.y,NumNecrosis)
            for i in range(0,len(xarray)):
                hasECMneighbor = False
                cell = self.cell_field[xarray[i], yarray[i], 0]
                if cell:
                    for neighbor, _ in self.get_cell_neighbor_data_list(cell):
                        if neighbor and neighbor.type == self.ECM: 
                            hasECMneighbor = True
                    if hasECMneighbor == False:
                        # What is point of this?
                        minDist = None
                        minCell = None
                        for cell in self.cell_list_by_type(self.ECM):
                            dist = np.sqrt((np.square(xarray[i]-cell.xCOM)) + (np.square(yarray[i] - cell.yCOM)))
                            if minDist == None or dist < minDist:
                                minDist = dist
                                minCell = cell
                        xarray[i] = cell.xCOM
                        yarray[i] = cell.yCOM
                        
                # Turn all fiber cells within certain radius into necrotic (as long as percentage doesn't exceed 30%)
                for cell in self.cell_list_by_type(self.FIBER):
                    r = np.sqrt((np.square(xarray[i]-cell.xCOM)) + (np.square(yarray[i] - cell.yCOM)))
                    if r < NecrosisRadius and NecroticPercentage < InitialNecroticPercentage:
                        cell.type = self.NECROTIC
                        NecroticVolume += cell.volume
                        NecroticPercentage = NecroticVolume / FiberVolume
        
        # Necrotic cells secrete HGF (proxy for nearby capillary burst bringing in HGF)
        for cell in self.cell_list_by_type(self.NECROTIC):
            HGFsecretor = self.get_field_secretor("HGF")
            HGFsecretor.secreteOutsideCellAtBoundary(cell, HGFsecretionVal)
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return

class NeutrophilSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # Initialize Neutrophils if starting with resident neutrophils
        
        # num_fiber_cells = len(self.cell_list_by_type(self.FIBER))
        # num_cells = np.round(num_fiber_cells * neuToInitialFiber)
        #print('Number of Neutrophils generated:')
        #print(num_cells)
        
        # for i in range(num_cells.astype(int)):
            # isPlaced = False
            
            # while not isPlaced:
                # xrand = np.random.randint(1,self.dim.x)
                # yrand = np.random.randint(1, self.dim.y)
                
                # cell = self.cellField[xrand,yrand,2]
                # if not cell:
                    # Neutrophil = self.new_cell(self.NEUTROPHIL)
                    # Neutrophil.targetVolume = 8.5
                    # Neutrophil.lambdaVolume = 50
                    # Neutrophil.dict['Birth'] = 0
                    # Neutrophil.dict['NumPhagocytosed'] = np.array([0])
                    # Neutrophil.dict['Apoptotic'] = 0
                    # Neutrophil.dict['secretionCountdown'] = -1
                    
                    # self.cell_field[xrand, yrand, 2] = Neutrophil
                    # isPlaced = True
        return

    def step(self, mcs):  
        # Recruit neutrophils if in correct time span
        if mcs < neutroStopRecruit: 
            if mcs > TimeNeuStart and (mcs - (TimeNeuStart)) % neutroRecruitmentFreq==0:   # make later so that it has a chance for MCP gradient to form from inital neutrophils 
                # Number of recuited neutrophils is proportional to number of necrotic cells
                fieldVal = [] 
                numRecruited = np.ceil(neutroRecruitmentProportion * len(self.cell_list_by_type(self.NECROTIC)))
                numRecruited = numRecruited.astype(int)
                
                # Find potential recruitment locations with low repel
                fieldRepel = self.field.REPEL
                coord = np.empty([0 , 2])
                numRecruitedTemp = numRecruited*recruitMultiplier 
                while numRecruitedTemp > 0:
                    xcoord1 = np.random.randint(1,self.dim.x)
                    ycoord1 = np.random.randint(1,self.dim.y)
                    coord1 = np.array([xcoord1, ycoord1])       
                    valueRepel = fieldRepel[xcoord1, ycoord1, 1] 
                    if valueRepel < 0.1:
                      coord = np.vstack((coord,coord1))
                      numRecruitedTemp -= 1               
                    
               # Prioritize locations based on hgf levels
                for x,y in coord:
                    fieldVal.append(self.field.HGF[x.item(), y.item(), 1])
                indexTopRegion = np.argsort(fieldVal)[-numRecruited:]  # sorts to get top numbers (going from the end because sorts low to high)
                toPlace = coord[indexTopRegion]
                
                # Place neutrophils if cell doesn't exist in that location
                for x,y in toPlace:
                    cell = self.cellField[x.item(),y.item(),1]
                    if not cell:
                        Neutrophil = self.new_cell(self.NEUTROPHIL)
                        self.cell_field[x.item(), y.item(), 1] = Neutrophil
                        Neutrophil.targetVolume = 8.5
                        Neutrophil.lambdaVolume = 50
                        #set neutrophil dictionary terms
                        Neutrophil.dict['Birth'] = mcs
                        Neutrophil.dict['NumPhagocytosed'] = 0
                        Neutrophil.dict['Apoptotic'] = 0
                        Neutrophil.dict['secretionCountdown'] = -1
                        #define neutrophil uptaking the hgf field
                        HGFsecretor = self.get_field_secretor("HGF") 
                        HGFsecretor.uptakeInsideCellAtBoundaryTotalCount(Neutrophil, neutrophilHGFUptakeMax, neutrophilHGFUptakeRel) #rel is a number between 0-1
              
        # Convert necrosis to low collagen ECM and secrete MCP and MMP when eating 
        MCPsecretor = self.get_field_secretor("MCP")
        MMPsecretor = self.get_field_secretor("MMP")
        if mcs > TimeNeuStart + NeutrophilPhagocytoseLagTime:        
            for cell in self.cell_list_by_type(self.NEUTROPHIL):
                CorrespondingNecrotic = self.cell_field[cell.xCOM, cell.yCOM, 0]
                # If cell below is necrotic convert to low collagen ecm
                if CorrespondingNecrotic and CorrespondingNecrotic.type == self.NECROTIC:
                    CorrespondingNecrotic.type = self.ECM
                    CorrespondingNecrotic.dict['collagen'] = 0.5
                    cell.dict['NumPhagocytosed'] += 1
                    MCPsecretor.secreteOutsideCellAtBoundary(cell, MCPsecretionVal)
                    MMPsecretor.secreteOutsideCellAtBoundary(cell, MMPsecretionVal)
                    cell.dict['secretionCountdown'] = necBreakdownLag # countdown to turn off mcp and mmp secretion
                    
        E2secretor = self.get_field_secretor("E2")
        for cell in self.cell_list_by_type(self.NEUTROPHIL):
            if cell.dict['Apoptotic'] == 0 and mcs%10 == 0:
                E2secretor.uptakeInsideCellAtBoundaryTotalCount(cell, 20/100.0, .3/100.0) # simulates E2 binding to ER on Neutrophils
            if cell.dict['secretionCountdown'] >= 0: # continue secreting
                MMPsecretor.secreteOutsideCellAtBoundary(cell, MMPsecretionVal)
                MCPsecretor.secreteOutsideCellAtBoundary(cell, MCPsecretionVal)
            else: # stop secretion
                MMPsecretor.secreteOutsideCellAtBoundary(cell, 0)
                MCPsecretor.secreteOutsideCellAtBoundary(cell, 0)
            cell.dict['secretionCountdown'] -= 1 # increment countdown

        #Make neutrophils apoptotic when appropriate 
        for cell in self.cell_list_by_type(self.NEUTROPHIL):
            if cell and (mcs - cell.dict['Birth']) >= neutroDeathTime or cell.dict['NumPhagocytosed'] >= 1:
                cell.dict['Apoptotic'] = 1 
                # stop movement
                cd = self.chemotaxisPlugin.getChemotaxisData(cell, "HGF")
                if cd:
                    cd.setLambda(0)

class SSCSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        # Track cell type or differentiation level of ssc for visualization
        self.track_cell_level_scalar_attribute(field_name='cellType', attribute_name='cellType')

    def start(self):
        # Initialize SSC (1 ssc per 4 fibers)
        num_fiber_cells = len(self.cell_list_by_type(self.FIBER))
        num_cells = np.round(num_fiber_cells * sscToInitialFiber)
        
        # Generate given number of cells where other cells don't already exist
        for i in range(num_cells.astype(int)):
            isPlaced = False
            
            while not isPlaced:
                xrand = np.random.randint(1,self.dim.x)
                yrand = np.random.randint(1, self.dim.y)
                
                cell = self.cellField[xrand,yrand,1]
                if not cell:
                    ssc = self.new_cell(self.SSC)
                    ssc.targetVolume = 8.5
                    ssc.lambdaVolume = 50
                    ssc.dict['activationState'] = 0 # 0 = nonactive, 1 = active
                    ssc.dict['time2activate'] = -1 # counter for activation lag
                    ssc.dict['time2divide'] = -1 # counter for division lag
                    ssc.dict['time2diff'] = -1 # counter for differentiation lag
                    ssc.dict['cellType'] = 0 # 0 = SSC, 1 = myoblast, 2 = myocyte
                    ssc.dict['moveCount'] = 0 # number of moves without being above low collagen ECM
                    
                    self.cell_field[xrand, yrand, 1] = ssc
                    isPlaced = True
                                  
                    # SSC chemotaxes toward MMP     
                    cd = self.chemotaxisPlugin.addChemotaxisData(ssc, "MMP")
                    cd.setLambda(3000.0)                   

    def step(self, mcs):
        # Calculate mean and max values of fields
        HGFfield = self.field.HGF
        MMPfield = self.field.MMP
        HGF_fieldval = np.zeros((xdim,ydim))
        MMP_fieldval = np.zeros((xdim,ydim))
        for i in range(0,self.dim.x-1):
            for j in range(0,self.dim.y-1):
                HGF_fieldval[i,j] = HGFfield[i,j,1]
                MMP_fieldval[i,j] = MMPfield[i,j,1]
        HGFMean = np.mean(HGF_fieldval)
        MMPMean = np.mean(MMP_fieldval)
        HGFmaxValue = HGFfield.max()  
        MMPmaxValue = MMPfield.max()
         
        # Recruit SSC proportional to sum of mean HGF and MMP every time step
        numRecruited = np.floor(recruitmentProportionSSC * (HGFMean + MMPMean))
        numRecruited = numRecruited.astype(int)

        # Randomly select spots and check to make sure potential spots have low enough repel
        fieldVal = []
        fieldRepel = self.field.REPEL
        coord = np.empty([0 , 2])
        numRecruitedTemp = numRecruited * recruitMultiplier
        while numRecruitedTemp > 0:
            xcoord1 = np.random.randint(1,self.dim.x)
            ycoord1 = np.random.randint(1,self.dim.y)
            coord1 = np.array([xcoord1, ycoord1])
            
            valueRepel = fieldRepel[xcoord1, ycoord1, 1] # Temp fix for layers
            if valueRepel < 0.1:
              coord = np.vstack((coord,coord1))
              numRecruitedTemp -= 1 
        
        # Prioritize locations with high HGF
        for x,y in coord:
            fieldVal.append(self.field.HGF[x.item(), y.item(), 1])
        indexTopRegion = np.argsort(fieldVal)[-numRecruited:]  # sorts to get top numbers (going from the end because sorts low to high)
        toPlace = coord[indexTopRegion]
        for x,y in toPlace:
            cell = self.cellField[x.item(),y.item(),1]
            if not cell:
                SSC = self.new_cell(self.SSC)
                self.cell_field[x.item(), y.item(), 1] = SSC
                SSC.targetVolume = 8.5
                SSC.lambdaVolume = 50
                SSC.dict['activationState'] = 1 # 0 = nonactive, 1 = active
                SSC.dict['time2activate'] = -1 # counter for activation lag
                SSC.dict['time2divide'] = -1
                SSC.dict['time2diff'] = -1
                SSC.dict['cellType'] = 0 # 0 = SSC, 1 = myoblast, 2 = myocyte
                SSC.dict['moveCount'] = 0
                
        # Actions
        E2secretor = self.get_field_secretor("E2")
        for cell in self.cell_list_by_type(self.SSC):
            # simulates E2 binding to ER on SSC
            E2secretor.uptakeInsideCellAtBoundaryTotalCount(cell, sscE2UptakeMax, sscE2UptakeRel)
            
            if cell.dict['activationState'] == 1 and cell.dict['time2divide'] < 0 and cell.dict['time2diff'] < 0 and cell.dict['time2activate'] < 0:
                
                cellBelow = self.cell_field[cell.xCOM, cell.yCOM, 0]
                                
                # no damage found in set number of moves -> deactivate
                if cellBelow and cellBelow.type != self.NECROTIC and cell.dict['moveCount'] > moveCountThreshold:
                    cell.dict['activationState'] = 0
                
                # if above low colagen ECM 
                if cellBelow and cellBelow.type == self.ECM and cellBelow.dict['collagen'] < lowCollagenCutoff:
                    # reset move counter
                    cell.dict['moveCount'] = 0 
                    
                    # check for neighboring fiber
                    fiberPresent = False
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cellBelow):
                        if neighbor and neighbor.type == self.FIBER:
                            fiberPresent = True
                            
                    # If not on fiber edge, migrate toward MMP
                    if not fiberPresent:
                        cd = self.chemotaxisPlugin.getChemotaxisData(cell, "MMP")
                        if cd:
                            lm = 3000.0
                            cd.setLambda(lm)
    
                    # If on fiber edge stop movement and divide or differentiate
                    if fiberPresent:    
                        # SSC or Myoblast dividing
                       if cell.dict['cellType'] == 1 or cell.dict['cellType'] == 0: 
                           # can make probabalistic
                            cell.dict['time2divide'] = sscDivisionTime
                        
                        # Myocte differentiating into fiber
                       if cell.dict['cellType'] == 2:
                            cell.dict['time2diff'] = sscDifferentiationTime
                        
                        # Stop moving while waiting to divide/diff    
                       cd = self.chemotaxisPlugin.getChemotaxisData(cell, "MMP")
                       if cd:
                            lm = 0.0
                            cd.setLambda(lm)
                            
                # if not above low ecm fiber move along mmp gradient
                else:
                    cd = self.chemotaxisPlugin.getChemotaxisData(cell, "MMP")
                    if cd:
                        lm = 3000.0
                        cd.setLambda(lm)
                    
                    cell.dict['moveCount'] += 1
                    
                    # apoptosis
                    if cell.dict['moveCount'] > sscApoptosisCutoff:
                        cell.targetVolume = 0
                        cell.lambdaVolume=1000000
                    
            elif cell.dict['activationState'] == 0:
                # attempt to activate by looking for growth factor
                field = self.field.HGF
                if field[cell.xCOM, cell.yCOM, 1] > sscActivationThreshold:
                    cell.dict['time2activate'] = sscActivationTime
                    
            # Increment lag counters 
            cell.dict['time2divide'] -= 1
            cell.dict['time2diff'] -=1  
            cell.dict['time2activate'] -=1

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return       

class SSCMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
    
    def step(self, mcs):
        # loop through all ssc and if lag times are over then divide or differentiate
        for cell in self.cell_list_by_type(self.SSC):
            # divide cell with possibility of changing type
            if cell.dict['time2divide'] == 0:      
                self.divide_cell_random_orientation(cell)
            
            # Turn low collage ecm into fiber
            if cell.dict['time2diff'] == 0: 
                cellBelow = self.cell_field[cell.xCOM, cell.yCOM, 0]
                if cellBelow and cellBelow.type == self.ECM and cellBelow.dict['collagen'] < lowCollagenCutoff:
                    cellBelow.type = self.FIBER
                    cell.targetVolume = 0   
            if cell.dict['time2activate'] == 0: 
                cell.dict['activationState'] = 1
    def update_attributes(self):
        # Clone cell
        self.clone_parent_2_child()   
        
        # 50 percent chance same type, 50 percent change next type (SSC->blast,blast->cyte)
        self.child_cell.dict['cellType'] = np.ceil(2*np.random.uniform()) -1 + self.parent_cell.dict['cellType']
        self.parent_cell.dict['cellType'] = np.ceil(2*np.random.uniform()) -1 + self.parent_cell.dict['cellType']
        
class MacrophageSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # Initialize resident macrophages
        num_fiber_cells = len(self.cell_list_by_type(self.FIBER))
        num_cells = np.floor((num_fiber_cells/256)*macToInitialFiber) # 3.7 from Kelley's code. Divide by 256 because 256 fiber elements in 1 fiber cell        
        
        for i in range(num_cells.astype(int)):
            isPlaced = False
            
            while not isPlaced:
                xrand = np.random.randint(1,self.dim.x)
                yrand = np.random.randint(1, self.dim.y) 
                # Check to make sure there is no cell in that random spot
                cell = self.cellField[xrand,yrand,1]
                if not cell:
                    Macrophage = self.new_cell(self.MACROPHAGE)
                    Macrophage.targetVolume = 8.5
                    Macrophage.lambdaVolume = 50
                    Macrophage.dict['Birth'] = 0
                    # MacrophageDeathTimeDistribution = np.random.randint(0.9*MacrophageDeathTime,1.1*MacrophageDeathTime)
                    MacrophagePoisson = -1* np.log(np.random.uniform()) *MacrophageDeathTime #poisson distribution
                    Macrophage.dict['DeathTime'] = 2000
                    Macrophage.dict['NumPhagocytosed'] = 0
                    Macrophage.dict['DeathIterator'] = MacrophageRandomDeath
                    Macrophage.dict['secretionCountdown'] = -1
                    self.cell_field[xrand, yrand, 1] = Macrophage
                    isPlaced = True
        return
        
    def step(self, mcs):
        # recruit macrophages if in the correct time interval
        if mcs < macroStopRecruit:
            if mcs > TimeMacroStart and (mcs - (TimeMacroStart)) % macroRecruitmentFreq == 0:   # make later so that it has a chance for MCP gradient to form from inital neutrophils 
                
                # Macrophage recruitment proportional to mean level of MCP
                MCPfield = self.field.MCP
                MCP_fieldval = np.zeros((xdim,ydim))
                for i in range(0,self.dim.x-1):
                    for j in range(0,self.dim.y-1):
                        MCP_fieldval[i,j] = MCPfield[i,j,1]
                
                MCPMean = np.mean(MCP_fieldval)
                numRecruited = np.ceil(origMacroRecruitmentProportion * MCPMean)
                numRecruited = numRecruited.astype(int)
                
                # Find grid squares with low repel
                fieldVal = []
                fieldRepel = self.field.REPEL
                coord = np.empty([0 , 2])
                numRecruitedTemp = numRecruited * recruitMultiplier
                while numRecruitedTemp > 0:
                    xcoord1 = np.random.randint(1,self.dim.x)
                    ycoord1 = np.random.randint(1,self.dim.y)
                    coord1 = np.array([xcoord1, ycoord1])
                    
                    valueRepel = fieldRepel[xcoord1, ycoord1, 1] 
                    if valueRepel < 0.1:
                      coord = np.vstack((coord,coord1))
                      numRecruitedTemp -= 1 
                
                # Prioritize locations with high HGF
                for x,y in coord:
                    fieldVal.append(self.field.HGF[x.item(), y.item(), 1])
                indexTopRegion = np.argsort(fieldVal)[-numRecruited:]  # sorts to get top numbers (going from the end because sorts low to high)
                toPlace = coord[indexTopRegion]
                for x,y in toPlace:
                    cell = self.cellField[x.item(),y.item(),1]
                    if not cell:
                        Macrophage = self.new_cell(self.MACROPHAGE)
                        self.cell_field[x.item(), y.item(), 1] = Macrophage
                        Macrophage.targetVolume = 8.5
                        Macrophage.lambdaVolume = 50
                        Macrophage.dict['Birth'] = mcs
                        Macrophage.dict['secretionCountdown'] = -1
                        # MacrophageDeathTimeDistribution = np.random.uniform(0.9*MacrophageDeathTime,1.1*MacrophageDeathTime)
                        MacrophagePoisson = -1* np.log(np.random.uniform()) *MacrophageDeathPoisson #poisson distribution
                        Macrophage.dict['DeathTime'] = MacrophagePoisson + MacrophageDeathTime
                        Macrophage.dict['NumPhagocytosed'] = 0
                        Macrophage.dict['DeathIterator'] = MacrophageRandomDeath
                        
        # After lag time allow macrophages to attempt to eat necrosis and neutrophils
        MMPsecretor = self.get_field_secretor("MMP")
        if mcs >= TimeMacroStart + MacrophagePhagocytoseLagTime:
            # turn necrotic cells below in low collagen ecm
            for cell in self.cell_list_by_type(self.MACROPHAGE): 
                CorrespondingNecrotic = self.cell_field[cell.xCOM, cell.yCOM, 0]
                if CorrespondingNecrotic and CorrespondingNecrotic.type == self.NECROTIC:
                    CorrespondingNecrotic.type = self.ECM
                    CorrespondingNecrotic.dict['collagen'] = 0.5
                    cell.dict['NumPhagocytosed'] += 1
                    MMPsecretor.secreteOutsideCellAtBoundary(cell, 100)
                    cell.dict['secretionCountdown'] = necBreakdownLag
            # Eat apoptotic neutrophils        
            for cell in self.cell_list_by_type(self.MACROPHAGE):
                for neighbor in self.get_cell_neighbor_data_list(cell):
                    if neighbor[0]:
                        if neighbor[0].type == self.NEUTROPHIL and neighbor[0].dict['Apoptotic'] == 1:
                            neighbor[0].targetVolume=0
                            neighbor[0].lambdaVolume=1000000  
                            cell.dict['NumPhagocytosed'] += 1    

        # define macrophage uptake of HGF, MCP, and E2
        HGFsecretor = self.get_field_secretor("HGF") 
        MCPsecretor = self.get_field_secretor("MCP")
        E2secretor = self.get_field_secretor("E2") 
        for cell in self.cell_list_by_type(self.MACROPHAGE):           
            HGFsecretor.uptakeInsideCellAtBoundaryTotalCount(cell, macrophageHGFUptakeMax, macrophageHGFUptakeRel) #rel is a number between 0-1
            MCPsecretor.uptakeInsideCellAtBoundaryTotalCount(cell, macrophageMCPUptakeMax, macrophageMCPUptakeMax)
            E2secretor.uptakeInsideCellAtBoundaryTotalCount(cell, macrophageE2UptakeMax, macrophageE2UptakeRel) # simulates E2 binding to ER on macrophages
        
        # If macrophage has lived its full lifespan then die 
        for cell in self.cell_list_by_type(self.MACROPHAGE):
            if mcs - cell.dict['Birth'] > cell.dict['DeathTime']:
                rand = random.randint(0,cell.dict['DeathIterator'])
                if rand == 0:
                    cell.targetVolume=0
                    cell.lambdaVolume=1000000
                if cell.dict['DeathIterator'] > 0:
                    cell.dict['DeathIterator'] -= 1
           
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return

class PlotSteppable(SteppableBasePy): # Create and add data to plots 
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # initialize graphs for CSA, #SSC, #Neutrophils, #Macro, and E2
        self.plot_win_csa = self.add_new_plot_window(title='Cross Sectional Area per step',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Area', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_win_csa.add_plot("FibArea", style='Lines', color='red', size=5)

        self.plot_win_ssc = self.add_new_plot_window(title='Stem Cell Count',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_win_ssc.add_plot("TotalSSC", style='Lines', color='red', size=5)
        self.plot_win_ssc.add_plot("ActiveSSC", style='Lines', color='orange', size=5)
        self.plot_win_ssc.add_plot("SSCNum", style='Lines', color='blue', size=5)
        self.plot_win_ssc.add_plot("MyoblastNum", style='Lines', color='purple', size=5)
        self.plot_win_ssc.add_plot("MyocyteNum", style='Lines', color='magenta', size=5)
        
        self.plot_win_neu = self.add_new_plot_window(title='Neutrophil Count',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_win_neu.add_plot("NeuNum", style='Lines', color='red', size=5)
        
        self.plot_win_macro = self.add_new_plot_window(title='Macrophage Count',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_win_macro.add_plot("MacroNum", style='Lines', color='red', size=5)
        
        self.plot_win_E2 = self.add_new_plot_window(title='Average E2',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Average E2', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_win_E2.add_plot("E2", style='Lines', color='red', size=5)
        return
        
    def step(self, mcs):
        # add data to graphs for CSA, #SSC, #Neutrophils, #Macro, and E2
        # arguments are (name of the data series, x, y)
        CellVolume = 0
        for cell in self.cell_list_by_type(self.FIBER):
            CellVolume += cell.volume
        self.plot_win_csa.add_data_point("FibArea", mcs, CellVolume)
              
        
        active = 0
        ssc = 0
        myoblast = 0
        myocyte = 0
        for cell in self.cell_list_by_type(self.SSC):
            if cell.dict['activationState'] == 1:
                active +=1
            if cell.dict['cellType'] == 0:
                ssc +=1
            if cell.dict['cellType'] == 1:
                myoblast +=1
            if cell.dict['cellType'] == 2:
                myocyte +=1
                  
        self.plot_win_ssc.add_data_point("TotalSSC", mcs, len(self.cell_list_by_type(self.SSC)))
        self.plot_win_ssc.add_data_point("ActiveSSC", mcs, active)
        self.plot_win_ssc.add_data_point("SSCNum", mcs, ssc)
        self.plot_win_ssc.add_data_point("MyoblastNum", mcs, myoblast)
        self.plot_win_ssc.add_data_point("MyocyteNum", mcs, myocyte)
     
        self.plot_win_neu.add_data_point("NeuNum", mcs, len(self.cell_list_by_type(self.NEUTROPHIL)))
           
        self.plot_win_macro.add_data_point("MacroNum", mcs, len(self.cell_list_by_type(self.MACROPHAGE)))
        
        # Calculate and plot E2 mean
        E2field = self.field.E2
        E2_fieldval = np.zeros((xdim,ydim))
        for i in range(0,self.dim.x-1):
            for j in range(0,self.dim.y-1):
                E2_fieldval[i,j] = E2field[i,j,0]

        E2Mean = np.mean(E2_fieldval)
        E2Mean = E2Mean.item()
        self.plot_win_E2.add_data_point("E2", mcs, E2Mean)
        
        # if mcs == 100:
            # if self.output_dir is not None:
                # output_path_CSA = Path(self.output_dir).joinpath('CSA.csv')
                # output_path_SSC = Path(self.output_dir).joinpath('StemCell.csv')
                # self.plot_win_csa.save_plot_as_data(output_path_CSA, CSV_FORMAT)
                # self.plot_win_ssc.save_plot_as_data(output_path_SSC, CSV_FORMAT)
                # # self.plot_win_csa.save_plot_as_data('CSA.csv')
                # # self.plot_win_csa.save_plot_as_data('StemCell.csv')

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        if self.output_dir is not None:
            output_path_CSA = Path(self.output_dir).joinpath('CSA.csv')
            output_path_SSC = Path(self.output_dir).joinpath('StemCell.csv')
            self.plot_win_csa.save_plot_as_data(output_path_CSA, CSV_FORMAT)
            self.plot_win_ssc.save_plot_as_data(output_path_SSC, CSV_FORMAT)
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return

class FileSteppable(SteppableBasePy): # Write output to files
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        # Define path to store data and open files with headers
        import os
        fileDir = os.path.dirname(os.path.abspath(__file__))
        fileDirUp = os.path.dirname(fileDir)
        
        fileName2 = fileDirUp+"\\logfile2MR.txt"
        self.file2 = open(fileName2,"a")
        outputText = "mcs, fiberVolume, sscCount, macroCount, neutroCount\n"
        self.file2.write(outputText)
        
        SSCfile = fileDirUp+"\\SSCdata.txt"
        self.file3 = open(SSCfile, "a")
        outputText = "mcs, cellID, cellType, xCOM, yCOM, activationState, cellTypeDifferentiation, HGF, MMP, MCP\n"        
        self.file3.write(outputText)
        
        neuFile = fileDirUp+"\\neuData.txt"
        self.file4 = open(neuFile, "a")
        outputText = "mcs, cellID, cellType, xCOM, yCOM, numPhagocytosed, Apoptotic, HGF, MMP, MCP\n"
        self.file4.write(outputText)
        
        macFile = fileDirUp+"\\macData.txt"
        self.file5 = open(macFile, "a")
        outputText = "mcs, cellID, cellType, xCOM, yCOM, numPhagocytosed, HGF, MMP, MCP\n"
        self.file5.write(outputText)
        
        CSAfile = fileDirUp+"\\CSAdata.txt"
        self.file6 = open(CSAfile, "a")
        outputText = "mcs, numHealthyFiber, numDamagedNecrotic, numRecoveringECM\n"
        self.file6.write(outputText)
        
        fieldFile = fileDirUp+"\\fieldData.txt"
        self.file7 = open(fieldFile, "a")
        outputText = "mcs, MCPMean, HGFMean, MMPMean\n"
        self.file7.write(outputText)
                
    def step(self, mcs):
        # Record fiber cell volume and agent counts
        CellVolume = 0
        for cell in self.cell_list_by_type(self.FIBER):
            CellVolume += cell.volume
        
        macroCount = len(self.cell_list_by_type(self.MACROPHAGE))
        sscCount = len(self.cell_list_by_type(self.SSC))
        neutroCount = len(self.cell_list_by_type(self.NEUTROPHIL))
        
        outputText = str(mcs)+", "+str(CellVolume)+", "+str(sscCount)+", "+str(macroCount)+", "+str(neutroCount)+"\n"
        self.file2.write(outputText)
        
        # Record fiber, ecm, necrosis counters
        numHealthy = len(self.cell_list_by_type(self.FIBER))
        numDamaged = len(self.cell_list_by_type(self.NECROTIC))
        global startingECM
        if mcs < 2:
            startingECM = len(self.cell_list_by_type(self.ECM))
        numRecovering = len(self.cell_list_by_type(self.ECM)) - startingECM 
        
        outputTextCSA = str(mcs)+", "+str(numHealthy)+", "+str(numDamaged)+", "+str(numRecovering)+"\n"
        self.file6.write(outputTextCSA)
        
        # Every 50 steps record ssc, neutrophil, macrophage, and field data
        if mcs %50 == 0: 
            for cell in self.cell_list_by_type(self.SSC):
                outputTextSSC = str(mcs)+", "+str(cell.id)+", "+str(cell.type)+", "+str(cell.xCOM)+", "+str(cell.yCOM)+ ", "+str(cell.dict['activationState'])+", " +str(cell.dict['cellType']) +", " +str(self.field.HGF[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MMP[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MCP[cell.xCOM, cell.yCOM,0])+"\n"
                self.file3.write(outputTextSSC)
                
            for cell in self.cell_list_by_type(self.NEUTROPHIL):
                outputTextNeu = str(mcs)+", "+str(cell.id)+", "+str(cell.type)+", "+str(cell.xCOM)+", "+str(cell.yCOM)+ ", "+str(cell.dict['NumPhagocytosed'])+", " +str(cell.dict['Apoptotic']) +", " +str(self.field.HGF[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MMP[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MCP[cell.xCOM, cell.yCOM,0])+"\n"
                self.file4.write(outputTextNeu)
                
            for cell in self.cell_list_by_type(self.MACROPHAGE):
                outputTextMac = str(mcs)+", "+str(cell.id)+", "+str(cell.type)+", "+str(cell.xCOM)+", "+ str(cell.yCOM) +", "+str(cell.dict['NumPhagocytosed'])+", " +str(self.field.HGF[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MMP[cell.xCOM, cell.yCOM,0]) +", " +str(self.field.MCP[cell.xCOM, cell.yCOM,0])+ "\n"
                self.file5.write(outputTextMac)
            
            MCPfield = self.field.MCP
            MCP_fieldval = np.zeros((xdim,ydim))
            for i in range(0,self.dim.x-1):
                for j in range(0,self.dim.y-1):
                    MCP_fieldval[i,j] = MCPfield[i,j,1]
            MCPMean = np.mean(MCP_fieldval)
            
            HGFfield = self.field.HGF
            HGF_fieldval = np.zeros((xdim,ydim))
            for i in range(0,self.dim.x-1):
                for j in range(0,self.dim.y-1):
                    HGF_fieldval[i,j] = HGFfield[i,j,1]
            HGFMean = np.mean(HGF_fieldval)
            
            MMPfield = self.field.MMP
            MMP_fieldval = np.zeros((xdim,ydim))
            for i in range(0,self.dim.x-1):
                for j in range(0,self.dim.y-1):
                    MMP_fieldval[i,j] = MMPfield[i,j,1]
            MMPMean = np.mean(MMP_fieldval)
            
            outputTextField = str(mcs)+", "+str(MCPMean)+", "+str(HGFMean)+", "+str(MMPMean)+"\n"
            self.file7.write(outputTextField)   
        
    def finish(self):
        # Display closing of files and close files
        self.file2.write("Closing file\n\n")
        self.file2.close()
        self.file3.write("Closing file\n\n")
        self.file3.close()
        self.file4.write("Closing file\n\n")
        self.file4.close()
        self.file5.write("Closing file\n\n")
        self.file5.close()
        self.file6.write("Closing file\n\n")
        self.file6.close()
        self.file7.write("Closing file\n\n")
        self.file7.close()
        return

    def on_stop(self):
        # Display closing of files and close files
        self.file2.write("Closing file\n\n")
        self.file2.close()
        self.file3.write("Closing file\n\n")
        self.file3.close()
        self.file4.write("Closing file\n\n")
        self.file4.close()
        self.file5.write("Closing file\n\n")
        self.file5.close()
        self.file6.write("Closing file\n\n")
        self.file6.close()
        self.file7.write("Closing file\n\n")
        self.file7.close()
        return  
  
  
class EstrogenSteppable(SteppableBasePy): # Cycle estrogen, secrete form capillaries, uptake by fibers
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):        
        return

    def step(self, mcs):
        # E2field = self.field.E2
        # E2_fieldval = np.zeros((self.dim.x,self.dim.y))
        # for i in range(0,self.dim.x-1):
            # for j in range(0,self.dim.y-1):
                # E2_fieldval[i,j] = E2field[i,j,0]
        # E2Tot = np.sum(E2_fieldval)    
        
        if estrogenOn == True:
            E2secretor = self.get_field_secretor("E2")        
            if mcs %9 == 0: # Every 6hrs is about 9 mcs if mcs = 41 mins
                global E2CycleIndex
                E2CycleIndex += 1 
                if E2CycleIndex >= len(cycleData):
                    E2CycleIndex = 0  
            
            systemicE2 = cycleData[E2CycleIndex] 
            
            # Secrete E2 from capillaries 
            for cell in self.cell_list_by_type(self.CAPILLARY):
                E2secretor.secreteOutsideCellAtBoundary(cell, systemicE2)                
                       
            # Estrogen is taken up by muscle cells
            for cell in self.cell_list_by_type(self.FIBER):
                # r = secretor.uptakeInsideCellTotalCount(cell, max, rel)
                E2secretor.uptakeInsideCellAtBoundaryTotalCount(cell, 20/100.0, .3/100.0) # simulates E2 binding to ER on fibers
      
     
     #Macrophage 
        # Estrogen promotes a less  severe infiltration of macrophages and improves macrophage phagocytosis (Liao et al., 2019)
        # global macroRecruitmentProportion
        # global MacrophagePhagocytoseLagTime
        # macroRecruitmentProportion = origMacroRecruitmentProportion * (MacroRecruitementProAlpha*E2Tot) # Decrease recruited proportion relative to E2
        # print("Total E2:")
        # print(E2Tot)
        # print("proportion:")
        # print(macroRecruitmentProportion)
        # MacrophagePhagocytoseLagTime -= (20*E2Tot) # Improved phagocytosis with E2 (shorter time to phagocytose)

        
        # #MCP Secretion
        # if mcs == 1:
            # global MCPsecretionVal
            # MCPsecretionVal *= (MCPAlpha/systemicE2)        
        
        #Neutrophil 
        # Estrogen increases infiltrating neutrophils 
        # global neutroRecruitmentProportion
        # neutroRecruitmentProportion *= (NeutroRecruitementProAlpha * (1/E2Mean))
        
        # #Necrotic
        # if mcs ==1:
            # global necBreakdownLag
            # necBreakdownLag *= (NecBreakdownLagAlpha * systemicE2)
        
        # #SSC
        # if mcs ==1:
            # global sscApoptosisCutoff
            # sscApoptosisCutoff *= (sscApoptosisAlpha * systemicE2)
        
        return

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return
        
class FiberSteppable(SteppableBasePy): # make fiber clusters and adjust repeal secretion
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        # Used to visualize clusters of fibers
        self.track_cell_level_scalar_attribute(field_name='clusterNum', attribute_name='clusterNum')

    def start(self):
        # join all connected fibers into one cluster
        fiberGroups = {}
        for cell in self.cell_list_by_type(self.FIBER): # Check all fibers
            cluster = None
            for neighbor, _ in self.get_cell_neighbor_data_list(cell): # Check current fiber neighbors
                if neighbor and neighbor.type == self.FIBER: # Make sure neighbor is a fiber
                    for group in fiberGroups: # For each cluster
                        if group != cluster:
                            if any(neighbor.id == item.id for item in fiberGroups[group]):
                            # if neighbor in fiberGroups[group]: Go through existing clusters to see if neighbor is in cluster
                                if cluster != None: # already assigned to a cluster 
                                    fiberGroups[cluster].extend(fiberGroups[group]) # Adding connecting groups
                                    fiberGroups.pop(group) # Gets rid of the group that we added so not duplicated
                                else: 
                                   fiberGroups[group].append(cell) # Not in a cluster so adding to nearby existing cluster
                                   cluster = group
                                break
            if cluster == None:
                fiberGroups[cell.clusterId] = [cell] # Starting a new cluster if no neighbors are a cluster
        
        compartmentNum = 1        
        for group in fiberGroups:
            for cell in fiberGroups[group]:
                self.reassign_cluster_id(cell, group)
                cell.dict['clusterNum'] = compartmentNum
                compartmentNum += 1
                
                # for compartments in self.clusters:
                    # for cell in compartments:
                        # print(cell.id)
                
        
 
       # join all connected fibers into one cluster
        # for cell in self.cell_list_by_type(self.FIBER):
            # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                # if neighbor and neighbor.type == self.FIBER:
                    # self.reassign_cluster_id(neighbor,cell.clusterId)
                    # for neighbor2, common_surface_area2 in self.get_cell_neighbor_data_list(neighbor):
                        # if neighbor2 and neighbor2.type == self.FIBER:
                            # self.reassign_cluster_id(neighbor2,cell.clusterId)
                            # for neighbor3, common_surface_area3 in self.get_cell_neighbor_data_list(neighbor2):
                                # if neighbor3 and neighbor3.type == self.FIBER:
                                    # self.reassign_cluster_id(neighbor3,cell.clusterId)
        
        return

    def step(self, mcs):
        #Make any fiber with one necrotic element permiable to cells
        secretor = self.get_field_secretor("REPEL")
        if mcs > 1: 
            # Loop through clusters
            for compartments in self.clusters:
                ifFiber = False
                ifNecrotic = False
                ifECM = False
                for cell in compartments:
                    if cell.type == self.FIBER:
                        ifFiber = True
                    if cell.type == self.NECROTIC:
                        ifNecrotic = True
                    if cell.type == self.ECM:
                        ifECM = True 
                
                # Make permeable
                # if ifFiber and ifNecrotic:
                    # for cell in compartments:
                        # secretor.secreteOutsideCellAtBoundary(cell, 0)
            
                # Make 'solid'
                if ifFiber and not ifNecrotic and not ifECM:
                    for cell in compartments:
                        secretor.secreteOutsideCellAtBoundary(cell, 1)
                    
class LymphaticSteppable(SteppableBasePy): # pulls cells toward lymph, removes chemokines, removes ssc, neutrophils, and macrophages
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        return

    def step(self, mcs):
        # Lymph uptakes hgf, mcp, mmp, e2, ssc, neutrophils, and macrophages
        HGFsecretor = self.get_field_secretor("HGF") #secretors for uptake
        MCPsecretor = self.get_field_secretor("MCP")
        MMPsecretor = self.get_field_secretor("MMP")
        E2secretor = self.get_field_secretor("E2")
        for lymp in self.cell_list_by_type(self.LYMPHATIC):
            # lymp uptakes fields
            HGFsecretor.uptakeInsideCellAtBoundaryTotalCount(lymp, 20, .8) #rel is a number between 0-1
            MCPsecretor.uptakeInsideCellAtBoundaryTotalCount(lymp, 20, .8)
            MMPsecretor.uptakeInsideCellAtBoundaryTotalCount(lymp, 20, .8)
            E2secretor.uptakeInsideCellAtBoundaryTotalCount(lymp, 20/100.0, .8/100.0)
            
            for neighbor, _ in self.get_cell_neighbor_data_list(lymp): # get all cells connected ot lymphatic cell
              if neighbor and (neighbor.type == self.SSC or neighbor.type == self.NEUTROPHIL or neighbor.type == self.MACROPHAGE):
                self.delete_cell(neighbor) #remove connected neighbors
            
            # force towards lymp, speed relative to distance from lymp
            for cell in self.cell_list_by_type(self.SSC, self.NEUTROPHIL, self.MACROPHAGE): 
                cell.lambdaVecX = lympForce/(cell.xCOM - lymp.xCOM) 
                cell.lambdaVecY = lympForce/(cell.yCOM - lymp.yCOM)
        return

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return


