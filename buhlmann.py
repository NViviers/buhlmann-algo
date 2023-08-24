import math

class Cmpt:
    HT: float
    M0: float
    DELTA_M: float

    def __init__(self, HT: float, M0: float, DELTA_M: float):
        self.HT = HT
        self.M0 = M0
        self.DELTA_M = DELTA_M

# Tissue compartments 1-16
COMPARTMENTS = [
    Cmpt(4, 106.4, 1.9082),
    Cmpt(8, 83.2, 1.5352),
    Cmpt(12.5, 73.8, 1.3847),
    Cmpt(18.5, 66.8, 1.278),
    Cmpt(27.0, 60.8, 1.2306),
    Cmpt(38.3, 55.6, 1.1857),
    Cmpt(54.3, 52.3, 1.1504),
    Cmpt(77.0, 50.1, 1.1223),
    Cmpt(109.0, 48.5, 1.0999),
    Cmpt(146.0, 47.2, 1.0844),
    Cmpt(187.0, 46.1, 1.0731),
    Cmpt(239.0, 45.1, 1.0635),
    Cmpt(305.0, 44.1, 1.0552),
    Cmpt(390.0, 43.1, 1.0478),
    Cmpt(498.0, 42.4, 1.0414),
    Cmpt(635.0, 41.8, 1.0359)
]

def cmpt_ppn2(compartment: int, ppn2: float, depth: int, time: float) -> float:
    cmpt = COMPARTMENTS[compartment]
    P = ppn2 + ((ppn2 * (depth / 10 + 1)) - ppn2) * (1 - 2 ** -(time / cmpt.HT))
    print(P)

cmpt_ppn2(0, 0.79, 40, 8.0)

# def cmpt_ppn2(compartment: int, ppn2: float, depth: int, time: float) -> float:
#     """Determines whether a decompression stop is required for a specific profile

#     Args:
#         compartment (int): Tissue compartment index
#         ppn2 (float): Partial pressure of nitrogen (in breathing gas) at surface
#         depth (float): Depth in meters
#         time (float): Bottom time at given depth

#     Returns:
#         bool: Whether a decompression stop is required
#     """
#     assert(compartment >= 0 and compartment <= 15)
#     cmpt = COMPARTMENTS[compartment]

#     P0 = (depth / 10 + 1) * ppn2
#     k = 0.693 / cmpt.HT

#     # Equilibrium pressure in fsw
#     A = cmpt.M0 + cmpt.DELTA_M * time

#     # PPN2 in tissue compartment in bar
#     P = ((P0 - A) * (1 - math.e ** (-k * time)) + A)
#     P_bar = P * 0.03048
#     print(P_bar)

#     return P

def cmpt_saturated(compartment: int, ppn2: float) -> bool:
    """Determines if a tissue compartment is saturated

    Args:
        compartment (int): Tissue compartment index
        ppn2 (float): Partial pressure of nitrogen in tissue compartment

    Returns:
        bool: Whether a compartment is satured or not
    """
    return ppn2 >= COMPARTMENTS[compartment].M0

# cmpt_ppn2(0, 0.79, 40, 7.0)
# cmpt_ppn2(0, 0.79, 40, 8.0)