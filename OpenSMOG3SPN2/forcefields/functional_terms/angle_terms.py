import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import sys
import os


def harmonic_angle_term(df_angles, use_pbc, force_group=2):
    angles = mm.HarmonicAngleForce()
    for i, row in df_angles.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        a3 = int(row['a3'])
        theta0 = row['theta0']
        k_angle = row['k_angle']
        angles.addAngle(a1, a2, a3, theta0, k_angle)
    angles.setUsesPeriodicBoundaryConditions(use_pbc)
    angles.setForceGroup(force_group)
    return angles






