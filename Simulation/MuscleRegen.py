
from cc3d import CompuCellSetup

from MuscleRegenSteppables import PassValuesSteppable
CompuCellSetup.register_steppable(steppable=PassValuesSteppable(frequency=1))

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
             
from MuscleRegenSteppables import MicrovesselSteppable
CompuCellSetup.register_steppable(steppable=MicrovesselSteppable(frequency=1))

from MuscleRegenSteppables import FibroblastSteppable
CompuCellSetup.register_steppable(steppable=FibroblastSteppable(frequency=1))

from MuscleRegenSteppables import FibroblastMitosisSteppable
CompuCellSetup.register_steppable(steppable=FibroblastMitosisSteppable(frequency=1))



        
from MuscleRegenSteppables import M2MitosisSteppable
CompuCellSetup.register_steppable(steppable=M2MitosisSteppable(frequency=1))

CompuCellSetup.run()
