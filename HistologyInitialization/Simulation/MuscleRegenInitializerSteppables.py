from cc3d.core.PySteppables import *
from cc3d.core.PySteppables import *
import numpy as np
from scipy.spatial import distance
import pandas as pd 

class MuscleRegenInitializerSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        fileName = 'add filepath to csv here'
        histoData = pd.read_csv(fileName, header=None)
        downscale = False
        if not downscale:
            for i in range(len(histoData)):
                x = histoData.iloc[i,0]
                y = histoData.iloc[i,1]
                tissueType = histoData.iloc[i,2]

                if tissueType == 3: # ECM 
                    cell_ecm = self.new_cell(self.ECM)
                    self.cellField[x.tolist(),y.tolist(),0] = cell_ecm

                if tissueType == 1: # ECM 
                    cell_ecm = self.new_cell(self.ECM)
                    self.cellField[x.tolist(),y.tolist(),0] = cell_ecm

        if downscale:
            for i in range(len(histoData)):
                x = histoData.iloc[i,0]
                y = histoData.iloc[i,1]
                tissueType = histoData.iloc[i,2]

                xAdjusted = int((x-1)/2 +1)
                yAdjusted = int((y-1)/2 +1)
                
           

                if tissueType == 3: # ECM
                    cell_ecm = self.new_cell(self.ECM)
                    for i in range(xAdjusted, xAdjusted+2):
                        for j in range(yAdjusted, yAdjusted+2):
                            self.cellField[i,j,0] = cell_ecm
                            cell_ecm.targetVolume = 10
                            cell_ecm.lambdaVolume = 1000

                if tissueType == 1: # ECM 
                    cell_ecm = self.new_cell(self.ECM)
                    for i in range(xAdjusted, xAdjusted+2):
                        for j in range(yAdjusted, yAdjusted+2):
                            self.cellField[i,j,0] = cell_ecm
                            cell_ecm.targetVolume = 10
                            cell_ecm.lambdaVolume = 1000

            # print(xAdjusted)
            # print(yAdjusted)
            # print('ecm volume', cell_ecm.volume)


    def step(self,mcs):
        if mcs == 1:
            foundSpace = False
            
            #place single seed to grow a cell
            while not foundSpace:
                xrand = np.random.uniform(0, self.dim.x) 
                yrand = np.random.uniform(0, self.dim.y)
               
                cell = self.cellField[xrand,yrand,0]
                if not cell:
                    cell_fib = self.new_cell(self.FIBER)
                    self.cellField[xrand, yrand,0] = cell_fib
                    cell_fib.targetVolume = 950000000
                    cell_fib.lambdaVolume = 1000
                    foundSpace = True
       
        allFilled = True
        # #if if there is empty space around a fiber cell (the seed from above), then find another random point and grow a fiber there too
        # for cell in self.cell_list_by_type(self.FIBER):
        #     for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
        #         if not neighbor:
        #             allFilled = False
                        
                
        if allFilled:
            foundSpace = False
            
            while not foundSpace:
                xrand = np.random.uniform(3, self.dim.x-3) 
                yrand = np.random.uniform(3, self.dim.y-3)
               
                cell = self.cellField[xrand,yrand,0]
                if not cell:
                    cell_fib = self.new_cell(self.FIBER)
                    self.cellField[xrand, yrand,0] = cell_fib
                    cell_fib.targetVolume = 950000000
                    cell_fib.lambdaVolume = 1000
                    foundSpace = True  

        
    def finish(self):
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return

threshold =  50
MinCapDist = 30
NumFibers = 1
class MakeCapillarySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        

    def start(self):
        # for x in range(0,self.dim.x):
        #     cell = self.cell_field[x, 0, 0]
        #     cell.type = self.WALL
            
        # for y in range(0,self.dim.y):
        #     cell = self.cell_field[0, y, 0]
        #     cell.type = self.WALL

        # for y in range(0, self.dim.y):
        #     cell = self.cell_field[358, y, 0]
        #     cell.type = self.WALL

        # for x in range(0, self.dim.x):
        #     cell = self.cell_field[x, 504, 0]
        #     cell.type = self.WALL


            
        return

    def step(self, mcs):
        if mcs == 3:
            CapillaryList = []
            for cell in self.cell_list_by_type(self.ECM):
                FibersNear = 0
                for fib in self.cell_list_by_type(self.FIBER):
                    dist = np.sqrt((np.square(fib.xCOM-cell.xCOM)) + (np.square(fib.yCOM - cell.yCOM)))
                    if dist < threshold:
                        FibersNear+=1
                
                NotNearby = True
                if FibersNear > NumFibers:
                    for cap in CapillaryList:
                        cap_dist = np.sqrt((np.square(cap.xCOM-cell.xCOM)) + (np.square(cap.yCOM - cell.yCOM)))
                        if cap_dist < MinCapDist:
                            NotNearby = False
                    if NotNearby:
                        CapillaryList.append(cell)
                        cell.type = self.CAPILLARY                    
                    
        return

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return

class MitosisSteppable(MitosisSteppableBase): # Only run to divide up fibers from orginial fiber segmentation image
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
    
    def start(self):
        # remove incomplete fibers
        for y in range(0,self.dim.y):
            cell = self.cellField[2.5,y,0]
            if cell is not None and cell.type == self.FIBER:
                cell.type = self.WALL

        for x in range(0,self.dim.x):
            cell = self.cellField[x,2.5,0]
            if cell is not None and cell.type == self.FIBER:
                cell.type = self.WALL

        for y in range(0,self.dim.y):
            cell = self.cellField[self.dim.x-1,y,0]
            if cell is not None and cell.type == self.FIBER:
                cell.type = self.WALL

        for x in range(0,self.dim.x):
            cell = self.cellField[x,self.dim.y-1,0]
            if cell is not None and cell.type == self.FIBER:
                cell.type = self.WALL
                    
        return
    
    def step(self, mcs):
        if mcs == 1 or mcs == 2:
            # remove incomplete fibers (done multiple times because of numerous compucell fiber seeds)
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL
            for cell in self.cell_list_by_type(self.WALL):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.WALL



        if mcs == 3:
            # Turn extended ECM into wall
            ecm = []
            wall = []
            wallCell = []
            for cell in self.cell_list_by_type(self.ECM):
                nextToFiber = False
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        nextToFiber = True
                if nextToFiber == False:
                    cell.type = self.WALL
                    wall.append((cell.xCOM, cell.yCOM))
                    wallCell.append(cell)
                else:
                    ecm.append((cell.xCOM, cell.yCOM))
                    
            localDist = 6
            dists = distance.cdist(wall, ecm).min(axis=1)
            for i in range(0,len(wallCell)):
                if dists[i] < localDist:
                    wallCell[i].type = self.ECM

        i = 0
        while i < 6: #randomly divide fiber fasicles into smaller regions for necrosis function
            cells_to_divide=[]
            for cell in self.cell_list_by_type(self.FIBER):
                if cell.volume > 80: # or 30 for smaller histology
                    cells_to_divide.append(cell)

            for cell in cells_to_divide:
                self.divide_cell_random_orientation(cell)
            i = i+1
             

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0 
        self.parent_cell.lambdaVolume/= 2.0
        self.parent_cell.fluctAmpl = 50
        self.clone_parent_2_child()        

        # for more control of what gets copied from parent to child use cloneAttributes function
        self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[])
        
class LymphaticSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        

    def start(self):

        return


    def step(self, mcs):
        if mcs == 4:
        
        # Location of lymphatic vessel 
            x = 25
            y = 170

            cell = self.cellField[x,y,0]
            
            cell.type = self.LYMPHATIC
            # print((len(self.cell_list_by_type(self.LYMPHATIC))))
            cell.targetVolume = 9500000
            cell.lambdaVolume = 1000
               
            for cell in self.cell_list_by_type(self.LYMPHATIC):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.FIBER:
                        neighbor.type = self.LYMPHATIC

   
        return 


    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return