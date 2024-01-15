from pymatgen.core import Lattice, Structure
import matgl


def to_kj_mol(ev_calculation):
    # convert ev/atom to kj/mol
    return round(ev_calculation * 96.4915666370759)


model = matgl.load_model("M3GNet-MP-2018.6.1-Eform")

gibbs_dict = {
    'enstatite': -349.394,
    'alpha-quartz': -204.646,
    'forsterite': -491.938,
    'fayalite': -329.668
}

"""
Enstatite
"""
struct = Structure.from_spacegroup("Pbca", Lattice.orthorhombic(18.2, 8.87, 5.20),
                                   ["Mg", "Mg", "Si", "Si", "O", "O", "O", "O", "O", "O"],
                                   [[.13, .33, .37],
                                    [.13, -.04, .37],
                                    [.03, -.35, .29],
                                    [.22, -.15, .04],
                                    [.06, .14, .20],
                                    [.06, .50, .20],
                                    [.05, -.25, .05],
                                    [.19, .35, .06],
                                    [.19, .01, .05],
                                    [.20, -.25, .30],
                                    ])
eform = model.predict_structure(struct)
print(f"The predicted formation energy for MgSiO3 is {to_kj_mol(float(eform.numpy()))} kj/mol.")
print(f"The formation energy for MgSiO3 from literature is {gibbs_dict['enstatite']} kj/mol.")


"""
Olivine (Forsterite)
"""

struct = Structure.from_spacegroup("Pbnm", Lattice.orthorhombic(4.7620, 10.2250, 5.9940),
                                   ["Mg", "Mg", "O", "O", "O", "Si"],
                                   [[0, 0, 0],
                                    [0.9896, 0.2776, 0.2500],
                                    [0.7667, 0.0918, 0.2500],
                                    [0.2202, 0.4477, 0.2500],
                                    [0.2781, 0.1633, 0.0337],
                                    [0.4226, 0.0945, 0.2500],
                                    ])
eform = model.predict_structure(struct)
print(f"The predicted formation energy for Mg2SiO4 is {to_kj_mol(float(eform.numpy()))} kj/mol.")
print(f"The formation energy for Mg2SiO4 from literature is {gibbs_dict['forsterite']} kj/mol.")

"""
Olivine (Fayalite)
"""

struct = Structure.from_spacegroup("Pbnm", Lattice.orthorhombic(4.79, 10.39, 6.06),
                                   ["Fe", "Fe", "O", "O", "O", "Si"],
                                   [[0, 0, 0],
                                    [0.9896, 0.2776, 0.2500],
                                    [0.7667, 0.0918, 0.2500],
                                    [0.2202, 0.4477, 0.2500],
                                    [0.2781, 0.1633, 0.0337],
                                    [0.4226, 0.0945, 0.2500],
                                    ])
eform = model.predict_structure(struct)
print(f"The predicted formation energy for Fe2SiO4 is {to_kj_mol(float(eform.numpy()))} kj/mol.")
print(f"The formation energy for Fe2SiO4 from literature is {gibbs_dict['forsterite']} kj/mol.")

"""
Alpha-Quartz
"""
struct = Structure.from_spacegroup(152, Lattice.trigonal(4.9137, 4.9137, 5.4047),
                                   ["Si", "O"],
                                   [[.4133, .2672, .1188],
                                    [.4697, 0.0, 0.0],
                                    ])
eform = model.predict_structure(struct)
print(f"The predicted formation energy for SiO2 is {to_kj_mol(float(eform.numpy()))} kj/mol.")
print(f"The formation energy for SiO2 from literature is {gibbs_dict['alpha-quartz']} kj/mol.")
