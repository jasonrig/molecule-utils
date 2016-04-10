class GamessJobFormatter(object):

    def __init__(self, molecule, basis_set="GBASIS=CCT", memory_replicated_gb=1, memory_distributed_gb=1):
        self.molecule = molecule
        self.basis_set = basis_set,
        self.mwords = memory_replicated_gb * 125
        self.memddi = memory_distributed_gb * 125