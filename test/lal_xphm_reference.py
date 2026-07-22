"""LALSuite v7.26.15 polarization oracle for the IMRPhenomXPHM Julia tests."""

import sys

import lal
import lalsimulation as lalsim


def main():
    values = [float(value) for value in sys.argv[1:]]
    if len(values) < 13:
        raise SystemExit(
            "expected m1 m2 s1x s1y s1z s2x s2y s2z distance_Gpc "
            "inclination phiRef fRef and at least one frequency"
        )

    m1, m2 = values[0:2]
    spins = values[2:8]
    distance_gpc, inclination, phi_ref, f_ref = values[8:12]
    frequencies_hz = values[12:]

    frequencies = lal.CreateREAL8Vector(len(frequencies_hz))
    frequencies.data[:] = frequencies_hz
    params = lal.CreateDict()
    lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion(params, 223)
    lalsim.SimInspiralWaveformParamsInsertPhenomXPConvention(params, 1)
    lalsim.SimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(params, 4)
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(params, 122022)
    # Disable both multibanding approximations in the numerical oracle.
    lalsim.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(params, 0.0)
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(params, 0.0)

    hp, hc = lalsim.SimIMRPhenomXPHMFrequencySequence(
        frequencies,
        m1 * lal.MSUN_SI,
        m2 * lal.MSUN_SI,
        *spins,
        distance_gpc * 1.0e9 * lal.PC_SI,
        inclination,
        phi_ref,
        f_ref,
        params,
    )
    for hp_i, hc_i in zip(hp.data.data, hc.data.data):
        print(
            f"{hp_i.real:.17e} {hp_i.imag:.17e} "
            f"{hc_i.real:.17e} {hc_i.imag:.17e}"
        )


if __name__ == "__main__":
    main()
