# https://en.wikipedia.org/wiki/B%C3%BChlmann_decompression_algorithm
import math

from typing import Optional

# Convenience class to ensure that we don't index incorrect values in a tuple or list for each compartment
class Cmpt:
    HT: float # Tissue half life
    M0: float # Saturation pressure
    DELTA_M: float # M-value gradient
    
	# The following are adjusted constants for the modified set C of Buhlmann's algorithm
    A: float # A constant for tissue compartment
    B: float # B constant for tissue compartment

    def __init__(self, HT: float, M0: float, DELTA_M: float, A: float, B: float):
        self.HT = HT
        self.M0 = M0
        self.DELTA_M = DELTA_M
        self.A = A
        self.B = B

# Tissue compartments 1-16
# M0 in fsw (feat sea water)
# Both fsw and bar should produce the same result
# These values are constants provided by ShearWater
# https://www.shearwater.com/wp-content/uploads/2019/05/understanding_m-values.pdf
COMPARTMENTS = [
    Cmpt(4, 106.4, 1.9082, 1.2599, 0.5050),
    Cmpt(8, 83.2, 1.5352, 1.0000, 0.6514),
    Cmpt(12.5, 73.8, 1.38478, 0.8618, 0.7222),
    Cmpt(18.5, 66.8, 1.278, 0.7562, 0.7825),
    Cmpt(27.0, 60.8, 1.2306, 0.6200, 0.8126),
    Cmpt(38.3, 55.6, 1.1857, 0.5043, 0.8434),
    Cmpt(54.3, 52.3, 1.1504, 0.4410, 0.8693),
    Cmpt(77.0, 50.1, 1.1223, 0.4000, 0.8910),
    Cmpt(109.0, 48.5, 1.0999, 0.3750, 0.9092),
    Cmpt(146.0, 47.2, 1.0844, 0.3500, 0.9222),
    Cmpt(187.0, 46.1, 1.0731, 0.3295, 0.9319),
    Cmpt(239.0, 45.1, 1.0635, 0.3065, 0.9403),
    Cmpt(305.0, 44.1, 1.0552, 0.2835, 0.9477),
    Cmpt(390.0, 43.1, 1.0478, 0.2610, 0.9544),
    Cmpt(498.0, 42.4, 1.0414, 0.2480, 0.9602),
    Cmpt(635.0, 41.8, 1.0359, 0.2327, 0.9653)
]

def cmpt_ppn2(compartment: int, fn2: float, depth: int, time: float) -> float:
	"""Computes the PPN2 of a tissue compartment with a given dive profile

	Args:
		compartment (int): Tissue compartment index (0-15)
		fn2 (float): FN2 of breathing gas
		depth (int): Depth of dive profile
		time (float): Bottom time of dive profile

	Returns:
		float: PPN2 of tissue compartment
	"""
	cmpt = COMPARTMENTS[compartment]

	return fn2 + ((fn2 * (depth / 10 + 1)) - fn2) * (1 - 2 ** -(time / cmpt.HT))

def cmpt_ceiling_pressure(compartment: int, ppn2: float) -> float:
	"""Computes the ceiling pressure of a tissue compartment

	Args:
		compartment (int): Tissue compartment index (0-15)
		ppn2 (float): PPN2 of tissue compartment

	Returns:
		float: Ceiling pressure of tissue compartment
	"""
	cmpt = COMPARTMENTS[compartment]
	
	return (ppn2 - cmpt.A) * cmpt.B

def cmpt_deco_needed(ppn2: float, ceiling_pressure: float) -> bool:
	"""Returns whether a tissue compartment's partial pressure exceeds it's ceiling pressure

	Args:
		ppn2 (float): PPN2 of tissue compartment
		ceiling_pressure (float): Ceiling pressure of tissue compartment

	Returns:
		bool: Whether a decompression stop is required
	"""
	return ppn2 >= ceiling_pressure

def cmpt_deco(compartment: int, ppn2: float, sf: int = 130) -> (float, float, float):
	"""Computes the time to decompress at the given depth

	Args:
		compartment (int): Tissue compartment index (0-15)
		ppn2 (float): Partial pressure of tissue compartment
		sf (int): Safety factor to scale the decompression by a percentage

	Returns:
		tuple: Duration of decompressing to do at 3 meters, 6 meters and 9 meters in minutes
	"""
	cmpt = COMPARTMENTS[compartment]
	k = 0.693 / cmpt.HT
	final = [0, 0, 0]

	# print(-1 * math.log((ppn2 - 1.3) / ppn2) / k)
	# after = 1.3 + (ppn2 - 1.3) * math.e ** (-k * 1.5)
	# print(f"Before: {ppn2}\nAfter: {after}")
	# print(-1 * math.log((ppn2 - 1.3) / after) / k)

	
	# Use a loop here instead, if 9m or 6m decom is less than one minute, then recalculate more time on the shallower depth
	final[2] = round(-1 * math.log((ppn2 - 1.9) / ppn2) / k) * (sf / 100)
	final[1] = round(-1 * math.log(0.3 / ppn2) / k) * (sf / 100)
	final[0] = round(-1 * math.log(0.3 / ppn2) / k) * (sf / 100)

	return tuple(final)

if __name__ == "__main__":
	tissue = 0
	fn2 = 0.79

	ppn2 = cmpt_ppn2(tissue, fn2, 40, 7)
	ceiling_pressure = cmpt_ceiling_pressure(tissue, ppn2)

	print(cmpt_deco(tissue, ppn2, ceiling_pressure))