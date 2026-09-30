"""
Identifying information for one completed CASCADE run.

    from cascade_pipeline.run_info import RunInfo

One object for the values every plotting function needs: name, folder, period, wave
height, sign convention. Details: scripts/cascade_pipeline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class RunInfo:
    """Identifying info for one completed CASCADE run.

    Attributes:
        run_name: Run name used in output filenames (e.g. run_name_hs).
        run_dir: Output directory for this run's files.
        start_year: Calendar year the run starts (period START_YEAR).
        end_year: Calendar year the run ends.
        Hs: Fixed significant wave height used for this run (m), shown in
            figure titles/legends. None omits it.
        flip_sign_model: CASCADE's x_s_TS increases landward (erosion).
            True (the convention used throughout cascade_pipeline's plotting)
            flips model output so a larger value means more seaward.
        background_erosion_on: Whether DOMAIN_BE_RATES has any non-zero
            entries for this run; used only for the "BE=on/off" label in
            figure titles and the GIF caption.
        model_name: Model name shown in figure titles/captions, e.g.
            "CASCADE". cascade_pipeline.shoreline is itself CASCADE-specific
            (it reads cascade.barrier3d / b3d.x_s_TS directly), so this
            defaults to "CASCADE" rather than a fully generic placeholder --
            override it if you're comparing against a different model run
            through the same figures.
        wave_climate: The run's wave settings as one line, e.g. "Hs 2.0 m,
            Tp 7.5 s, asym 0.6, high-angle 0.5". Shown in the subtitle only
            when RateComparisonConfig.show_wave_climate is set. None omits it.
    """

    run_name: str
    run_dir: str
    start_year: int
    end_year: int
    Hs: float = None
    flip_sign_model: bool = True
    background_erosion_on: bool = True
    model_name: str = "CASCADE"
    wave_climate: str = None
