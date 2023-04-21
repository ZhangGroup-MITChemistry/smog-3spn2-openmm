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


def dna_3spn2_stacking_term(df_stackings, use_pbc, force_group=7):
    stackings = mm.CustomCompoundBondForce(3, f"""energy;
                energy=rep+f2*attr;
                rep=epsilon*(1-exp(-alpha*(dr)))^2*step(-dr);
                attr=epsilon*(1-exp(-alpha*(dr)))^2*step(dr)-epsilon;
                dr=distance(p2,p3)-sigma;
                f2=max(f*pair2,pair1);
                pair1=step(dt+{np.pi}/2)*step({np.pi}/2-dt);
                pair2=step(dt+{np.pi})*step({np.pi}-dt);
                f=1-cos(dt)^2;
                dt=rng*(angle(p1,p2,p3)-theta0);""")
    stackings.addPerBondParameter('epsilon')
    stackings.addPerBondParameter('sigma')
    stackings.addPerBondParameter('theta0')
    stackings.addPerBondParameter('alpha')
    stackings.addPerBondParameter('rng')
    # add parameters
    for i, row in df_stackings.iterrows():
        parameters = [row['epsilon'], row['sigma'], row['theta0'], row['alpha'], row['rng']]
        a1, a2, a3 = int(row['a1']), int(row['a2']), int(row['a3'])
        stackings.addBond([a1, a2, a3], parameters)
    stackings.setUsesPeriodicBoundaryConditions(use_pbc)
    stackings.setForceGroup(force_group)
    return stackings



