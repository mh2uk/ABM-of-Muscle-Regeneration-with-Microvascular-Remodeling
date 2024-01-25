
from cc3d import CompuCellSetup

 
# start with just this uncommented 
from MuscleRegenInitializerSteppables import MuscleRegenInitializerSteppable
CompuCellSetup.register_steppable(steppable=MuscleRegenInitializerSteppable(frequency=1))
  
# Once above is run and piff is saved and added to xml, comment the above and uncomment the below steppables      

# from MuscleRegenInitializerSteppables import MitosisSteppable
# CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

# from MuscleRegenInitializerSteppables import MakeCapillarySteppable
# CompuCellSetup.register_steppable(steppable=MakeCapillarySteppable(frequency=1))
        
# from MuscleRegenInitializerSteppables import LymphaticSteppable
# CompuCellSetup.register_steppable(steppable=LymphaticSteppable(frequency=1))

CompuCellSetup.run()
