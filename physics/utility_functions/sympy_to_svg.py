import sympy as smp


def sympyToSvg(smpFunction, filename: str):
    return smp.preview(smpFunction, output="svg", viewer="file", filename=filename)
