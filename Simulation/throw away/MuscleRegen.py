
from cc3d import CompuCellSetup
        

# from MuscleRegenSteppables import ImportInitializerSteppable
# CompuCellSetup.register_steppable(steppable=ImportInitializerSteppable(frequency=1))

from MuscleRegenSteppables import ConstraintInitializerSteppable
CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))

# from MuscleRegenSteppables import MitosisSteppable
# CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

from MuscleRegenSteppables import FiberSteppable
CompuCellSetup.register_steppable(steppable=FiberSteppable(frequency=1))

from MuscleRegenSteppables import PlotSteppable
CompuCellSetup.register_steppable(steppable=PlotSteppable(frequency=1))

from MuscleRegenSteppables import FileSteppable
CompuCellSetup.register_steppable(steppable=FileSteppable(frequency=1))

from MuscleRegenSteppables import EstrogenSteppable
CompuCellSetup.register_steppable(steppable=EstrogenSteppable(frequency=1))
        
from MuscleRegenSteppables import NecrosisSteppable
CompuCellSetup.register_steppable(steppable=NecrosisSteppable(frequency=1))

        
from MuscleRegenSteppables import NeutrophilSteppable
CompuCellSetup.register_steppable(steppable=NeutrophilSteppable(frequency=1))


from MuscleRegenSteppables import SSCSteppable
CompuCellSetup.register_steppable(steppable=SSCSteppable(frequency=1))


from MuscleRegenSteppables import MacrophageSteppable
CompuCellSetup.register_steppable(steppable=MacrophageSteppable(frequency=1))

from MuscleRegenSteppables import SSCMitosisSteppable
CompuCellSetup.register_steppable(steppable=SSCMitosisSteppable(frequency=1))
       
from MuscleRegenSteppables import LymphaticSteppable
CompuCellSetup.register_steppable(steppable=LymphaticSteppable(frequency=1))

CompuCellSetup.run()
