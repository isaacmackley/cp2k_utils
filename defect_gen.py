from pymatgen.analysis.defects.generators import ChargeInterstitialGenerator
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure
from pymatgen.io.common import VolumetricData
from pymatgen.io.cif import CifWriter



def int_gen(cube_file, specie, output_filename, int_no=10):
    """
    """

    volumetric_data = VolumetricData.from_cube(cube_file)
    chgcar = Chgcar(volumetric_data.structure, volumetric_data.data)

    ints = [[],[]]

    for defect in ChargeInterstitialGenerator(max_insertions=int_no).generate(chgcar, insert_species=specie):
        ints[0].append(defect.site.specie)
        ints[1].append(defect.site.frac_coords)
    

    structures = {}

    i=0
    while i < len(ints[0]):
        structure = volumetric_data.structure
        structures[i+1] = structure.append(species=ints[0][i],coords=ints[1][i])
        CifWriter(structures[i+1]).write_file(f'{output_filename}{i+1}.cif')
        i=i+1